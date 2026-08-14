# create a multiple sequence alignment and segment the protein 
# into conserved and less conserved sections
# Author: Klaus Foerstemann

# install libraries


if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("Biostrings")) BiocManager::install("Biostrings")  
if (!require("msa")) BiocManager::install("msa")
if (!require("ggtree")) BiocManager::install("ggtree") 
if (!require("seqinr")) {install.packages("seqinr")} 
if (!require("Peptides")) {install.packages("Peptides")} 
if (!require("treeio")) {install.packages("treeio")} 
if (!require("ape")) {install.packages("ape")} 
if (!require("tidyverse")) {install.packages("tidyverse")}

library("tidyverse")
library("ggtree")
library("msa")
library("ape")
library("treeio")
library("Biostrings")
library("seqinr")
library("Peptides")

# define the size of the sliding window up here for convenience
# large window sizes lead to fewer sections, small window sizes give more detail/resolution
# the ideal compromise will depend on the particular alignment and the goal of the analysis

window_size <-30   # suggested default is 30; can go as high as hundred if appropriate; a minimum of 13 is enforced
# values below 13 will throw an error downstream (division by 25 & rounded for threshold, must not be 0 though)
if (window_size<13) {window_size <- 13}


# define parameters required later

aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # aminoacids



# read in the multi-line fasta file
file_list <- list.files(pattern = "_unique_BLAST_results\\.fa$")



for(fasta_filename in file_list) {
  # remove the extension
  fasta_filename_noext <- sub("(\\..*)", "", fasta_filename)
  # remove _unique_BLAST_results if present
  fasta_filename_noext <- strsplit(fasta_filename_noext, "_unique_BLAST_results")[[1]][1]
  
  #create a new subdirectory to place the analysis files
  starting_point <- format(Sys.time(), "%y%m%d_%H-%M")
  subdir_name <- paste(fasta_filename_noext, starting_point, sep="_")
  dir.create(subdir_name, showWarnings = FALSE)
  
  multi_fasta_sequences <- readAAStringSet(file=fasta_filename)
  

  
  
  ###################################################
  # generate MSA, segment and extract the sequences #
  ###################################################
  
  
  # generate the alignment; algorithm options are ClustalW, ClustalOmega and Muscle
  # see the documentation of the msa package for algorithm details and additional parameters
  # https://ftp.gwdg.de/pub/misc/bioconductor/packages/3.17/bioc/vignettes/msa/inst/doc/msa.pdf
  #
  # You may need to adjust the gap opening and extension penalties!
  filename_alignment <- paste0(subdir_name,"/msa_", fasta_filename_noext, ".fasta")
  
  print("calculating MSA and distance matrix, this can take a short moment.")
  
  generated_MSA <- msa(multi_fasta_sequences, 
                       "Muscle",        # Muscle is recommended because it understands ambiguous AA codes
                       type="protein",
                       gapOpening=30,   # set this rather high to avoid too many gaps
                       gapExtension=0.5
  )
  
  # build a phylogenetic tree of the sequences
  # Convert to seqinr alignment for distance calculation
  MSA_seqinr <- msaConvert(generated_MSA, type = "seqinr::alignment")
  # Compute the distance matrix
  MSA_dist_matrix <- dist.alignment(MSA_seqinr, "similarity")
  
  
  
  # Build a Neighbor-Joining tree
  sequence_tree <- bionjs(MSA_dist_matrix)
  
  # Plot the tree circular
  sequence_tree_plot <- ggtree(sequence_tree, options(ignore.negative.edge=TRUE), layout="circular")+
    geom_tiplab(aes(angle=angle), size=240/length(sequence_tree$tip.label)) + 
    theme_tree() +
    ggtitle(label= paste("MSA of ", fasta_filename_noext)) +
    theme(plot.title = element_text(hjust = 0.5))   # centered position
  # show the plot
  print(sequence_tree_plot) 
  # save the plot
  filename <- paste0(subdir_name,"/sequence_tree_circ_", fasta_filename_noext, ".pdf")
  ggsave(filename, sequence_tree_plot, width=8, height=10)
  
  # Plot the tree linear
  sequence_tree_plot <- ggtree(sequence_tree, options(ignore.negative.edge=TRUE), layout="rectangular")+
    geom_tiplab(aes(angle=0), size=240/length(sequence_tree$tip.label)) + 
    theme_tree() +
    ggtitle(label= paste("MSA of ", fasta_filename_noext)) +
    theme(plot.title = element_text(hjust = 0.5))   # centered position
  # show the plot
  print(sequence_tree_plot) 
  # save the plot
  filename <- paste0(subdir_name,"/sequence_tree_lin_", fasta_filename_noext, ".pdf")
  ggsave(filename, sequence_tree_plot, width=8, height=10)
  
  # extract the consensus sequence and abstract it to 0 and 1 (1 = conserved residue)
  consensus_sequence_msa <- msaConsensusSequence(generated_MSA)
  abstracted_consensus <- gsub("[aA-zZ]", "1", gsub("[-?]", "0", consensus_sequence_msa))
  # convert this to a vector of numbers
  abstracted_consensus_numeric <- as.numeric(unlist(strsplit(abstracted_consensus, "")))
  
  # The resulting positions need to be combined into segments for further analysis. 
  # A sliding window approach works reasonably well here, more sophisticated procedures can certainly be implemented.
  
  # calculate a sliding window defined above, sum the abstracted values, round to 0 or 1 to smooth the curve
  
  thresholding <- round(window_size/25, 0) # choose the parameters window size above so that you get a reasonable segmentation
  smooth_consensus <- abstracted_consensus_numeric # initialize new vector
  # ramp up until the first window size
  for (ramp_position in 1:window_size) {
    smooth_consensus[ramp_position] <- round(sum(abstracted_consensus_numeric[seq(1:ramp_position)])/thresholding, 0)
  }
  # then continue until the end
  for (sliding_window_position in (window_size + 1):length(abstracted_consensus_numeric)) {
    window_start <- sliding_window_position - window_size
    if(round(sum(abstracted_consensus_numeric[window_start:sliding_window_position])/thresholding, 0)>=1) smooth_consensus[(window_start+thresholding):sliding_window_position] <- 1
    else smooth_consensus[(sliding_window_position-thresholding):sliding_window_position] <- 0
  }
  
  
  
  
  
  # convert the profile into sections of conserved (1) and non-conserved (0) sequence
  sequence_sections <- vector()
  current_state <- smooth_consensus[1]
  first_pos <- 1
  last_pos <- 1
  section_counter <- 1
  polished_consensus <- vector() # create a copy for polishing
  for (i in 1: length(smooth_consensus)) {
    if (smooth_consensus[i]==current_state) {
        next
      }
    # this part will be executed when a change in the smoothed sequence occurs
    last_pos <- i-1 # marks the end of the section in the smoothed sequence;
    print(paste("rough pass last:", last_pos))
    # now look at this in detail and take the state change from the non-smoothed sequence
    fine_pass_last <- last_pos # initialize with the rough value
    # polish the end of the section
    start_pos <- 1 #initialize with 1
    # Analysis starts at position 1 of the sequence or further downstream, as log as that is beyond the window size
    if (last_pos > window_size) {start_pos = last_pos - window_size} #change only if this does not lead to a negative number
    for (polish_region in (start_pos:(last_pos + window_size))) {
      
      if (abstracted_consensus_numeric[polish_region]==current_state) {next}
      # once a change occurs, this part is executed
      print(paste("polishing section change position:", polish_region))
      fine_pass_last <- polish_region -1 # set the beginning of the new section
      print(paste("fine pass last:", fine_pass_last))
      break # stop going further in the polishing loop (not elegant but works)
    }
    
    # now make sure the section is at least 2 aminoacids long
    if ((fine_pass_last - first_pos) >=2) {
      sequence_sections[section_counter] <- paste0(first_pos,":",fine_pass_last)
      names(sequence_sections)[section_counter] <- current_state
      polished_consensus[first_pos:fine_pass_last] <- current_state
      print(paste("old state:", current_state))
      current_state <- smooth_consensus[i] # set new state according to the state change detected in the smoothed sequence
      print(paste("new state:", current_state))
      section_counter <- section_counter +1
      first_pos <- fine_pass_last+1
      print(paste("fine pass first adjusted", first_pos))
    }
  }
  
  # finish off the last segment
  last_pos <- length(smooth_consensus)
  sequence_sections[section_counter] <- paste0(first_pos,":",last_pos)
  polished_consensus[first_pos:last_pos] <- current_state
  names(sequence_sections)[section_counter] <- current_state
  
  # generate a plot to document the segmentation
  abstracted_data_df <- data.frame(c(seq(1:length(abstracted_consensus_numeric))), abstracted_consensus_numeric, polished_consensus)
  
  abstracted_consensus_plot <- ggplot(abstracted_data_df, aes(x = abstracted_data_df[,1], y=abstracted_data_df[,2])) +
    geom_point(shape=1) +
    geom_line(y = abstracted_data_df[,3], color = "green", linewidth = 1) +
    theme_classic() +
    theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=14, hjust = 0.5), axis.text.x = element_text(size=18), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
    labs(title = paste("Segments in", fasta_filename_noext), x = "MSA position", y= "0=variable, 1=conserved") 
  
  # show the plot
  print(abstracted_consensus_plot)
  # save the plot
  filename <- paste0(subdir_name,"/segmentation_", fasta_filename_noext, ".pdf")
  ggsave(filename, abstracted_consensus_plot, width=12, height=8)
  
  
  
  
  # now use this to name the positions in the MSA for each sequence.
  sequence_set <- as.character(unmasked(generated_MSA))
  ID_index <- vector()
  accession_index <- vector()
  for(i in 1: length(sequence_set)) {
    ID_index <- c(ID_index, (i*2)-1 )
    accession_index <- c(accession_index, i*2)
  }
  
  sequence_names <- unlist(str_split(rownames(generated_MSA), "\\|"))[accession_index] # pipe symbol needs escape for matching
  sequence_species_ID <- unlist(str_split(rownames(generated_MSA), "\\|"))[ID_index]
  
  # generate names for the dataframe columns
  state_dictionary <- c("variable", "conserved") # just need to add 1 to the section values to get the corresponding word
  DF_column_names <- c("species_ID", "accession")
  for (i in 1:length(sequence_sections)) {
    print(state_dictionary[as.numeric(names(sequence_sections[i]))+1])
    DF_column_names <- c(DF_column_names, paste(i, state_dictionary[as.numeric(names(sequence_sections[i]))+1], sep="_"))
  }
  sectioned_sequences <- as.data.frame((matrix("", ncol= length(sequence_sections)+2, nrow= length(sequence_names))))
  colnames(sectioned_sequences) <- DF_column_names
  sectioned_sequences$accession <- sequence_names
  sectioned_sequences$species_ID <- sequence_species_ID
  
  # generate numeric versions of the sequence sections
  section_starts <- vector()
  section_ends <- vector()
  
  for (i in 1:length(sequence_sections)) {
    numbers <- strsplit(sequence_sections[i], ":")
    section_starts[i] <- as.numeric(numbers[[1]][1])
    section_ends[i] <- as.numeric(numbers[[1]][2])
  }
  
  # go through every sequence
  for (i in 1:length(sequence_set)) {
    # split the sequence into individual positions
    split_sequence <- strsplit(sequence_set[i], "")
    # and go through every section
    for (j in 1:length(sequence_sections)) {
      sectioned_sequences[i,j+2] <- paste(split_sequence[[1]][section_starts[j]:section_ends[j]], collapse="")
    }
  }
  
  # copy into new dataframe, then remove all dashes
  
  collapsed_sectioned_sequences <- sectioned_sequences
  for (i in 2: ncol(collapsed_sectioned_sequences)) {
    collapsed_sectioned_sequences[,i] <-  gsub("-", "", collapsed_sectioned_sequences[,i])
  }
  
  # save the extracted sequence segments for future use
  
  filename_gapped <- paste0(subdir_name, "/gapped_", fasta_filename_noext, ".csv")
  write.csv(sectioned_sequences, file=filename_gapped)
  
  filename_collapsed <- paste0(subdir_name, "/collapsed_", fasta_filename_noext, ".csv")
  write.csv(collapsed_sectioned_sequences, file=filename_collapsed)
  
  
  
  ###################################
  # analysis of sequence properties #
  ###################################
  
  # compute size distribution for each fragment
  size <- vector()
  section <- vector()
  type <- vector()
  section_length_df <- data.frame(size, section, type)
  section_length_table <- collapsed_sectioned_sequences # copy this to replace the sequences by their lengths
  for (i in 3:length(collapsed_sectioned_sequences)) {
    for (j in 1:length(collapsed_sectioned_sequences[,i])) {
      section_length_df <- rbind(section_length_df, c(nchar(collapsed_sectioned_sequences[j,i]), i, as.character(colnames(collapsed_sectioned_sequences[i])))) # for plot
      section_length_table[j,i] <-as.numeric(nchar(collapsed_sectioned_sequences[j,i]) )  # to save in table
    }
  }
  colnames(section_length_df) <- c("size", "section", "type")
  filename_size_df <- paste0(subdir_name, "/sizes_", fasta_filename_noext,".csv")
  write.csv(section_length_table, file=filename_size_df)
  
  sizes_plot <- ggplot(section_length_df, aes(y=as.numeric(section_length_df[,1]), x=as_factor(section_length_df[,3]))) +   # note that as_factor() leaves the order untouched, while as.factor() will reorder alphabetically; not desired for >10 sections
    labs(title = paste("Sizes in", fasta_filename_noext), subtitle = paste(length(collapsed_sectioned_sequences[,1]), " sequences analyzed"), x = "sequence section", y= "length(AA)") +
    theme_classic() +
    theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=14, hjust = 0.5), axis.text.x = element_text(size=18, angle = 90, hjust=0, vjust = 0.5), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
    
    geom_jitter(shape=1, alpha=1, color="darkgrey") + 
    geom_boxplot(width=0.7, linewidth=0.5, color="black", fill="salmon", alpha=0.5, outliers = FALSE) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y),0)), vjust = -1, color = "black", size = 4)
  
  # show the plot
  print(sizes_plot)
  # save the plot
  filename <- paste0(subdir_name,"/size_analysis_", fasta_filename_noext, ".pdf")
  ggsave(filename, sizes_plot, width=10, height=8)
  
  # compute GRAVY score (hydrophobicity) for each fragment
  GRAVY <- vector()
  section <- vector()
  type <- vector()
  section_gravy_df <- data.frame(size, section, type)
  section_gravy_table <- collapsed_sectioned_sequences # copy this to replace the sequences by their lengths
  for (i in 3:length(collapsed_sectioned_sequences)) {
    for (j in 1:length(collapsed_sectioned_sequences[,i])) {
      current_gravy <- hydrophobicity(seq = collapsed_sectioned_sequences[j,i], scale = "KyteDoolittle")
      section_gravy_df <- rbind(section_gravy_df, c(current_gravy, i, as.character(colnames(collapsed_sectioned_sequences[i])))) # for plot
      section_gravy_table[j,i] <-as.numeric(current_gravy)  # to save in table
    }
  }
  colnames(section_gravy_df) <- c("gravy", "section", "type")
  filename_size_df <- paste0(subdir_name, "/gravy_", fasta_filename_noext,".csv")
  write.csv(section_gravy_table, file=filename_size_df)
  
  gravy_plot <- ggplot(section_gravy_df, aes(y=as.numeric(section_gravy_df[,1]), x=as_factor(section_gravy_df[,3]))) +
    labs(title = paste("GRAVY in", fasta_filename_noext), subtitle = paste(length(collapsed_sectioned_sequences[,1]), " sequences analyzed"), x = "sequence section", y= "GRAVY score (Kyte-Doolittle)") +
    theme_classic() +
    theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=14, hjust = 0.5), axis.text.x = element_text(size=18, angle = 90, hjust=0, vjust = 0.5), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
    
    geom_jitter(shape=1, alpha=1, color="darkgrey") + 
    geom_boxplot(width=0.7, linewidth=0.5, color="black", fill="salmon", alpha=0.5, outliers = FALSE) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)), vjust = -1, color = "black", size = 4)
  
  
  # show the plot
  print(gravy_plot)
  # save the plot
  filename <- paste0(subdir_name,"/gravy_analysis_", fasta_filename_noext, ".pdf")
  ggsave(filename, gravy_plot, width=10, height=8)
  
  
  # compute the pI of each fragment
  pI_analysis <- data.frame(cbind(collapsed_sectioned_sequences$species_ID, collapsed_sectioned_sequences$accession))
  for (i in 3:ncol(collapsed_sectioned_sequences)){
    #print(i)
    for (j in 1:length(collapsed_sectioned_sequences[,i])) {
      sequence_pI <- s2c(collapsed_sectioned_sequences[j,i]) # break up the sequence string into individual AA
      if (length(sequence_pI)>=1) {pI_analysis[j,i] <- computePI(sequence_pI)}
      
    }
  }
  colnames(pI_analysis) <- DF_column_names
  # save this table
  filename_pI_df <- paste0(subdir_name, "/isoelectric_", fasta_filename_noext,".csv")
  write.csv(pI_analysis, file=filename_pI_df)
  
  
  # combine into new df for plotting
  isoelectric <- vector()
  section <- vector()
  type <- vector()
  pI_sections <- data.frame(isoelectric, section, type)
  for (i in 3:length(pI_analysis)) {
    for (j in 1:length(pI_analysis[,i])) {
      pI_sections <- rbind(pI_sections, c(as.numeric(pI_analysis[j,i]), i, as.character(colnames(collapsed_sectioned_sequences[i]))))
      
    }
  }
  colnames(pI_sections) <- c("isoelectric", "section", "type")
  # make the plot
  pI_plot <- ggplot(pI_sections, aes(y=as.numeric(pI_sections[,1]), x=as_factor(pI_sections[,3]))) +
    labs(title = paste("pI in",sub("(\\..*)", "", fasta_filename)), subtitle = paste(length(collapsed_sectioned_sequences[,1]), " sequences analyzed"), x = "sequence section", y= "isoelectric point") +
    theme_classic() +
    theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=14, hjust = 0.5), axis.text.x = element_text(size=18, angle = 90, hjust=0, vjust = 0.5), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
    
    geom_jitter(shape=1, alpha=1, color="darkgrey") + 
    geom_boxplot(width=0.7, linewidth=0.5, color="black", fill="salmon", alpha=0.5, outliers = FALSE) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)),vjust = -1, color = "black", size = 4)
  
  
  # show the plot
  print(pI_plot)
  # save the plot
  filename <- paste0(subdir_name,"/pI_analysis_", fasta_filename_noext, ".pdf")
  ggsave(filename, pI_plot, width=10, height=8)
  
  
  
  # calculate AA statistics
  AA_analysis_percent <- list("accession" = collapsed_sectioned_sequences$accession)
  AA_analysis_percent[[2]] <- list("species_ID" = collapsed_sectioned_sequences$species_ID)
  for (i in 3:ncol(collapsed_sectioned_sequences)){
    #print(i)
    current_section <- list()
    names(current_section[1]) <- colnames(collapsed_sectioned_sequences[,i])
    for (j in 1:length(collapsed_sectioned_sequences[,i])) {
      sequence_AA <- s2c(collapsed_sectioned_sequences[j,i])
      
      if (length(sequence_AA)>=1) {current_section[j] <- list(round(100*seqinr::count(sequence_AA, wordsize=1, alphabet=aa, freq=TRUE),2))}
      else {current_section[j] <- NA }
    }
    AA_analysis_percent[[i]] <- current_section
  }
  names(AA_analysis_percent) <- DF_column_names
  
 
  # create one matrix per section and fill that with the AA percentages
  # get the column names from the first set that actually has a sequence
  n <- 1
  exist_names <- 0
  while(exist_names<2) {
    n <- n+1
    exist_names <- length(AA_analysis_percent[3][[1]][[n]])
    if (n>600) {break} # include an emergency stop just in case
  } 
  
  composition_colnames <- names(AA_analysis_percent[3][[1]][[n]])
  for (i in 3:length(AA_analysis_percent)) {
    composition_df <- data.frame(cbind(unlist(AA_analysis_percent$species_ID), (unlist(AA_analysis_percent$accession)), matrix(ncol=length(composition_colnames))))
    rownames(composition_df) <- NULL
    
    for (j in 1:length(AA_analysis_percent[[i]])) {
      for (k in 1:length(composition_colnames)) {
        composition_df[j,k+2] <- as.numeric(AA_analysis_percent[[i]][j][[1]][k])
      }
    }
    colnames(composition_df) <- c("accession", "species_ID", composition_colnames)
    
    # save each table separately
    filename_composition_df <- paste0(subdir_name, "/composition_section_", (i-2), "_", fasta_filename_noext,".csv")
    write.csv(composition_df, file=filename_composition_df) 
    
    # create a plot for this section
    composition_pivot <- pivot_longer(composition_df, cols=colnames(composition_df[3:length(composition_df)]), names_to = "Aminoacid")
    AA_plot <- ggplot(composition_pivot, aes(y=as.numeric(composition_pivot$value), x=composition_pivot$Aminoacid)) +
      labs(title = paste("AA frequency ", fasta_filename_noext), subtitle = paste("Sequence section", (i-2), "present in", length(composition_df[which(!is.na(composition_df[,3])),3]), "sequences"), x= "", y= "percent") +
      theme_classic() +
      theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=18, hjust = 0.5), axis.text.x = element_text(size=18, angle = 90, hjust=0, vjust = 0.5), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
      expand_limits(y=0) +
      geom_jitter(shape=1, alpha=1, color="darkgrey") + 
      geom_boxplot(width=0.7, linewidth=0.5, color="black", fill="salmon", alpha=0.5, outliers = FALSE) 
    
    # show the plot
    print(AA_plot)
    # save the plot
    filename <- paste0(subdir_name,"/AA_comp_sec_", (i-2),"_", fasta_filename_noext, ".pdf")
    ggsave(filename, AA_plot, width=10, height=8)
    
  }
  
  
  # calculate "interesting" dipeptide frequencies
  
   
  
  dipeptide_freq_selected <- list("accession" = collapsed_sectioned_sequences$accession)
  dipeptide_freq_selected[[2]] <- list("species_ID" = collapsed_sectioned_sequences$species_ID)
  for (i in 3:ncol(collapsed_sectioned_sequences)){
    current_section <- list()
    names(current_section[1]) <- colnames(collapsed_sectioned_sequences[,i])
    for (j in 1:length(collapsed_sectioned_sequences[,i])) {
      sequence_AA <- s2c(collapsed_sectioned_sequences[j,i])
      if (length(sequence_AA)>1) {current_section[j] <- list(round(100*seqinr::count(sequence_AA, wordsize=2, alphabet=aa, freq=TRUE),2))}
      else {current_section[j] <- NA }
    }
    dipeptide_freq_selected[[i]] <- current_section
  }
  names(dipeptide_freq_selected) <- DF_column_names
  
  #
  
  # create one matrix per section and fill that with the dipeptide frequencies
  
  dipeptide_colnames <- names(dipeptide_freq_selected[3][[1]][[n]])
  for (i in 3:length(dipeptide_freq_selected)) {
    dipeptide_df <- data.frame(cbind(unlist(dipeptide_freq_selected$species_ID), (unlist(dipeptide_freq_selected$accession))), matrix(ncol = length(dipeptide_colnames)))
    rownames(dipeptide_df) <- NULL
    
    for (j in 1:length(dipeptide_freq_selected[[i]])) {
      for (k in 1:length(dipeptide_colnames)) {
        if(length(dipeptide_freq_selected[[i]][j][[1]]) > 1) {dipeptide_df[j,k+2] <- as.numeric(dipeptide_freq_selected[[i]][j][[1]][k])}
      }
    }
    
    
    colnames(dipeptide_df) <- c("accession", "species_ID", dipeptide_colnames)
    
    # save each table separately
    filename_dipeptide_df <- paste0(subdir_name, "/dipeptides_section_", (i-2), "_", fasta_filename_noext,".csv")
    write.csv(dipeptide_df, file=filename_dipeptide_df) 
    
    # calculate the average for each dipeptide
    dipep_average <- vector()
    for (k in 1:length(dipeptide_colnames)) {
      dipep_average[k] <- mean(as.numeric(unlist(dipeptide_df[k+2])), na.rm=TRUE)
    }
    # get the top 40 dipeptides
    dipep_top_40 <- vector()
    # get the top dipeptide values, maximum 30 but only those that are non-zero
    dipep_top_40 <- order(as.numeric(dipep_average), decreasing = TRUE)[1:min(length(dipep_average[which(dipep_average>0)]),40)]
    dipep_df_columns <- c(1,2,(dipep_top_40+2)) # get the name and ID column as well
    #dipeptide_df_selected <- data.frame()
    dipeptide_df_selected <- data.frame(dipeptide_df[dipep_df_columns])
    
    # create a plot for this section
    # define what to plot
    
    
    # comment out one of these lines to choose the sorting of the graph:
    
    # invoke this line to sort the top 40 dipeptides alphabetically (aminoacids)
    #dipeptide_pivot <- pivot_longer(dipeptide_df_selected, cols=colnames(dipeptide_df_selected[3:length(dipeptide_df_selected)]), names_to = "Dipeptides", values_drop_na = TRUE)
    
    # invoke this line to sort the dipeptides according to the median in decreasing order (i.e. the most frequent to the left)
    dipeptide_pivot <- pivot_longer(dipeptide_df_selected, cols = -c(accession, species_ID), names_to = "Dipeptides", values_drop_na = TRUE) %>% mutate(Dipeptides = fct_reorder(Dipeptides, value, .fun = median, .desc = TRUE ))
   
    # now generate the plot
     dipeptide_plot <- ggplot(dipeptide_pivot, aes(y=as.numeric(dipeptide_pivot$value), x=dipeptide_pivot$Dipeptides)) +
      labs(title = paste("Top 40 Dipeptide frequency ", fasta_filename_noext), subtitle = paste("Sequence section", (i-2), "present in", length(dipeptide_df_selected[which(!is.na(dipeptide_df_selected[,3])),3]), "sequences"), x= "", y= "percent") +
      theme_classic() +
      theme(plot.title = element_text(size = 24, hjust = 0.5), plot.subtitle = element_text(size=18, hjust = 0.5), axis.text.x = element_text(size=18, angle = 90, hjust=0, vjust = 0.5), axis.text.y=element_text(size=18), axis.title.y = element_text(size=24), axis.title.x = element_text(size=24)) +
      expand_limits(y=0) +
      geom_jitter(shape=1, alpha=1, color="darkgrey") + 
      geom_boxplot(width=0.7, linewidth=0.5, color="black", fill="salmon", alpha=0.5, outliers = FALSE) 
    
    # show the plot
    print(dipeptide_plot)
    # save the plot
    filename <- paste0(subdir_name,"/Dipep_comp_sec_", (i-2),"_", fasta_filename_noext, ".pdf")
    ggsave(filename, dipeptide_plot, width=10, height=8)
    
  }
}
