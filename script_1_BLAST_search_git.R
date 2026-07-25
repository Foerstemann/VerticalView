#### BLAST search for homologous proteins based on a given protein template ####
#### Author: Klaus Foerstemann


#########################################################################################
# This script accesses the NCBI BLAST server. PLEASE RESPECT THEIR GUIDELINES!          #
# That will also help to make your requests run faster - they will cut you off rapidly! #
# If you embark on a higher volume project, a local BLAST server may be preferable.     #
#########################################################################################

# From https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#developerinfo 
#
# ...
# Usage Guidelines
#
# The NCBI BLAST servers are a shared resource. We give priority to interactive users. 
# In order to ensure availability of the service to the entire community, we may limit searches for some high volume users. 
# Interactive users of the NCBI webpages through a web browser should not encounter problems.
# We will move searches of users who submit more than 100 searches in a 24 hour period to a slower queue, or, in extreme cases, will block the requests. 
# To avoid problems, API users should comply with the following guidelines:
#  
# Do not contact the server more often than once every 10 seconds.
#
# Do not poll for any single RID more often than once a minute.
#
# Use the URL parameter email and tool, so that the NCBI can contact you if there is a problem.
#
# Run scripts weekends or between 9 pm and 5 am Eastern time on weekdays if more than 50 searches will be submitted.
# ...
#


# install and load libraries

if (!require("curl")) install.packages("curl")
if (!require("rvest")) install.packages("rvest")
if (!require("rentrez")) install.packages("rentrez")
if (!require("taxize")) {install.packages("taxize")}
if (!require("protr")) {install.packages("protr")} 
if (!require("treeio")) {install.packages("treeio")} 
if (!require("ape")) {install.packages("ape")} 
if (!require("ggplot2")) install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require("ggtree")) BiocManager::install("ggtree") 

library("ggtree")
library("protr") 
library("taxize") 
library("curl")
library("rvest")
library("rentrez")
library("treeio")
library("ape")
library("ggplot2")


##########################
# set desired parameters #
##########################


# Define the entrez API key - this is not strictly required but recommended.
# You can create this in your pubmed account,
# then store it e.g. in the environment function pubmed_key().
# Of course you can also directly code it here as a string: set_entrez_key("YOUR_KEY"),
# but then make sure to remove it before sharing the script with others.
# Simply comment out the line if you are not using an entrez API key.
  set_entrez_key(pubmed_key())


# define the query protein sequence(s)
# this is a named list and may contain multiple proteins if desired
protein_query_list <- list(Dm_Ago2 = "MGKKDKNKKGGQDSAAAPQPQQQQKQQQQRQQQPQQLQQPQQLQQPQQLQQPQQQQQQQPHQQQQQSSRQQPSTSSGGSRASGFQQGGQQQKSQDAEGWTAQKKQGKQQVQGWTKQGQQGGHQQGRQGQDGGYQQRPPGQQQGGHQQGRQGQEGGYQQRPPGQQQGGHQQGRQGQEGGYQQRPSGQQQGGHQQGRQGQEGGYQQRPPGQQQGGHQQGRQGQEGGYQQRPSGQQQGGHQQGRQGQEGGYQQRPPGQQQGGHQQGRQGQEGGYQQRPPGQQQGGHEQGRQGQEGGYQQRPSGQQQGGHQQGRQGQEGGYQQRPSGQQQGGHQQGRQGQEGGYQQRPSGQQQGGHQQGRQGQEGGYQQRPPGQQPNQTQSQGQYQSRGPPQQQQAAPLPLPPQPAGSIKRGTIGKPGQVGINYLDLDLSKMPSVAYHYDVKIMPERPKKFYRQAFEQFRVDQLGGAVLAYDGKASCYSVDKLPLNSQNPEVTVTDRNGRTLRYTIEIKETGDSTIDLKSLTTYMNDRIFDKPMRAMQCVEVVLASPCHNKAIRVGRSFFKMSDPNNRHELDDGYEALVGLYQAFMLGDRPFLNVDISHKSFPISMPMIEYLERFSLKAKINNTTNLDYSRRFLEPFLRGINVVYTPPQSFQSAPRVYRVNGLSRAPASSETFEHDGKKVTIASYFHSRNYPLKFPQLHCLNVGSSIKSILLPIELCSIEEGQALNRKDGATQVANMIKYAATSTNVRKRKIMNLLQYFQHNLDPTISRFGIRIANDFIVVSTRVLSPPQVEYHSKRFTMVKNGSWRMDGMKFLEPKPKAHKCAVLYCDPRSGRKMNYTQLNDFGNLIISQGKAVNISLDSDVTYRPFTDDERSLDTIFADLKRSQHDLAIVIIPQFRISYDTIKQKAELQHGILTQCIKQFTVERKCNNQTIGNILLKINSKLNGINHKIKDDPRLPMMKNTMYIGADVTHPSPDQREIPSVVGVAASHDPYGASYNMQYRLQRGALEEIEDMFSITLEHLRVYKEYRNAYPDHIIYYRDGVSDGQFPKIKNEELRCIKQACDKVGCKPKICCVIVVKRHHTRFFPSGDVTTSNKFNNVDPGTVVDRTIVHPNEMQFFMVSHQAIQGTAKPTRYNVIENTGNLDIDLLQQLTYNLCHMFPRCNRSVSYPAPAYLAHLVAARGRVYLTGTNRFLDLKKEYAKRTIVPEFMKKNPMYFV",
                           Dm_Nup54 = "MSFFGSNTSLGATSTPAKTTGGLFGSPFGGTAATSQPAPAFGAQATSTPAFGAQPATSAFGAGSAFGATAAAPAFGAATGTSAFGGSAFGSTPAFGAATTTTAGTGLGGGGFGGFGAAPATSQAGLFGAPATSAAPPAFSGFGQQAAASTAPASGFSGFGTTTTSAPAFGGFGTSQSTGFGGGAFGSTFGKPANTTVTPGFGGFGGTSFMLGQPQQQPAPISADEAFAQSILNVSIFGDERDKIVAKWNYLQATWGTGKMFYSQSAAPVDITPENVMCRFKAIGYSRMPGKDNKLGLVALNFCRELSAVKPHQQQVIQTLHSLFGSKPNMLVHIDSIKELENKKCQIVIYVEEKLQHAPNESKRILATELSNYLNQATLKPQLNNLGVVEALALVLPDEDQLREYLENPPRGVDPRMWRQANSDNPDPTLYIPVPMVGFNDLKWRVKCQEQETDTHALYIKKVESELTELKKRHATATAKILEHKRKLAELSHRILRIIVKQECTRKVGTSLTPEEEALRTKLQNMLAVVSAPTQFKGRLSELLSQMRMQRNQFAANGGAEYALDKEAEDEMKTFLTMQQRAMEVLSDTVNKDLRALDVIIKGLPELRQS",
                           Dm_SC35 = "MSNGGGAGGLGAARPPPRIDGMVSLKVDNLTYRTTPEDLRRVFERCGEVGDIYIPRDRYTRESRGFAFVRFYDKRDAEDALEAMDGRMLDGRELRVQMARYGRPSSPTRSSSGRRGGGGGGGSGGRRRSRSRSPMRRRSRSPRRRSYSRSRSPGSHSPERRSKFSRSPVRGDSRNGIGSGSGGLAPAASRSRSRS"
)

# keep the longest protein isoform per species or the first one in the list?
keep_longest = FALSE   # set this to TRUE if desired


# maximum number of retrieved BLAST hit sequences (this is for the non-unique set!)
max_hits <- 500    # a minumum of 20 is suggested because the phylogenetic analysis needs at least 3 distinct species in the unique dataset; a reasonable value for a deep analysis is 500; going beyond that may lead to server issues   

# define the desired BLAST parameters - this corresponds pretty much to the interactive NCBI BLAST webpage
# first try the parameters out interactively on the NCBI Blast webpage rather than here in the script
# once a good parameter set has been found, plug the values in the script here and let it run

BLAST_db <- "refseq_protein"               # "refseq_protein" is the only one that will work for the taxonomy analysis
BLAST_filter <- "T"                        # low-complexity filter; default for BLAST is "F" i.e. disabled. For proteins with IDR's: Use "T" or "L" to enable; option "mF" or "mT" or "mL" will influence at lookup selectively
BLAST_expect <- 0.01                       # default suggested is 0.01 - use lower for more stringent and higher for less stringent search; beware of CPU overtime if you increase this beyond 0.1 !
BLAST_word <- 5                            # default suggested is 5, allowed values are only 3, 5 and 6; a value of 3 can lead to CPU overtime errors
BLAST_gap <- paste("11", "2", sep="%20")   # existence and extension cost of gaps; default is "11 2", allowed values are "11 2", "10 2" , "9 2", "8 2", "7 2", "6 2", "13 1", "12 1", "11 1", "10 1", "9 1"  
                                           # the space between the two numbers must be encoded by %20 in the URL 
BLAST_matrix <- "BLOSUM62"                 # default is "BLOSUM62", other options are "BLOSUM45", "BLOSUM50", "BLOSUM80", "BLOSUM90", "PAM250","PAM30" and "PAM70"
BLAST_composition <- 2                     # algorithm to use for composition based statistics; default is 2, allowed values are 0, 1, 2, and 3. 
BLAST_short <- "false"                     # automatically adjust parameters for short query sequences; default is "false" but may be set to "true" if desired.
BLAST_alignments <- 1                      # we will not use the alignments, so keep this low to save data transfer




################
# main routine #
################


# loop over the query list given above
for (query_counter in 1:length(protein_query_list)) {
 
  
  
  query_seq <- gsub("\\s", "", protein_query_list[[query_counter]]) # remove all whitespace characters from the query sequence
  
  # BLAST search

  # construct the URL with the query sequence
  url_request <- paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastp&FORMAT_TYPE=TEXT",
                        "&DATABASE=", BLAST_db,
                        "&FILTER=", BLAST_filter,
                        "&EXPECT=", BLAST_expect,
                        "&WORD_SIZE=", BLAST_word,
                        "&GAPCOSTS=", BLAST_gap,
                        "&MATRIX=", BLAST_matrix,
                        "&COMPOSITION_BASED_STATISTICS=", BLAST_composition,
                        "&SHORT_QUERY_ADJUST=", BLAST_short,
                        "&HITLIST_SIZE=", max_hits,
                        "&ALIGNMENTS=", BLAST_alignments,
                        "&QUERY=",query_seq)
  
  # submit the request 
 

  response <- curl_fetch_memory(url_request)
  content <- read_html(response$content)
  # extract the ID
  request_id_element <- html_node(content, "#rid")
  request_id <- html_attr(request_id_element, "value") 
  print(paste("NCBI BLAST request ID:", request_id))
  
  # construct the URL with the query sequence
 
  url_results <- paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=TEXT&ALIGNMENT_VIEW=Tabular",
                        "&DESCRIPTIONS=", max_hits,
                        "&ALIGNMENTS=", BLAST_alignments,
                        "&RID=", request_id)
  
  # poll for results every 60 seconds by checking whether the parsing of the response worked
  BLAST_hits <- vector()
  while (length(BLAST_hits)<1){
    
    results_response <- curl_fetch_memory(url_results)
    content_response <- read_html(results_response$content)
    response_parsed <- html_text(content_response)
    response_split <- unlist(strsplit(response_parsed,"\n"))
    if (length(response_split[grep("Time since submission", response_split)])) {
       print(response_split[grep("Time since submission", response_split)])
     } else if (length(response_split[grep("CPU usage limit", response_split)])) {
       print(response_parsed)
       break
     }
      
    BLAST_hits <- grep("Query_", response_split)
    # wait 60 seconds between status checks
    system.time({ Sys.sleep(60) })
  }

  accession_list <- vector()
  for (line in BLAST_hits){
    #print(line)
    new_accession <- strsplit(response_split[line], "\t")[[1]][2]
    accession_list <- c(accession_list, new_accession) 
  }
  
  protein_sequences <- vector()
  
  for (record_counter in 1:length(accession_list)) {
    # assign this directly instead of appending to the list via c(...) saves memory
    protein_sequences[record_counter] <- entrez_fetch(db = "protein", id = accession_list[record_counter], rettype = "fasta")
    system.time({ Sys.sleep(0.1) }) # introduce a minimal waiting time to not overload the server; improves stability
  }
  
  # save this dataset immediately as a fasta file so that it can be used later  
  # it can take a while to reach this point and one should be able to pick up from here if necessary.
  proteins_to_save <- noquote(protein_sequences)
  proteins_to_save2 <- gsub(">", "\n>" , proteins_to_save) # this is necessary to remove a whitespace character at the beginning of the line
  filename_BLAST_results <- paste0(names(protein_query_list[query_counter]),"_BLAST_results.fa")
  # using sink is unusual but necessary to prevent the quotes from appearing in the output file
  sink(file=filename_BLAST_results, append=FALSE, type="output", split = FALSE)
  cat(proteins_to_save2)
  sink(file = NULL)
  
  
  # to continue the work: pretend this is a file and parse it into a list with the command from the protr package
  protein_fasta_seqs <- readFASTA(textConnection(proteins_to_save2))
  
  
  # extract the species name between the square brackets
  species_names <- sub(".*\\[(.*?)\\].*", "\\1", protein_sequences)
  
  # construct a short version thereof
  # initialize the new list
  short_spec_names <- vector(mode = "list", length = length(species_names))
  # split the species names on the space
  split_species_names <- strsplit(species_names, " ")
  
  
  for (i in 1:length(species_names)) {
    # cobine the first letter of the genus with the first three letters of the species
    short_spec_names[[i]] <- paste0(substr(split_species_names[[i]][1], 1,1), substr(split_species_names[[i]][2], 1, 3))
  }
  
  # combine everything into a dataframe for convenience
  annotated_BLAST_results <- data.frame(critters=species_names, accession=paste(short_spec_names, names(protein_fasta_seqs), sep="_"), sequence=trimws(unlist(unname(protein_fasta_seqs))))
  
  # now keep only the one sequence for each critter;
  # by default this will be the first sequence in the list for each species
  # but you can un-comment the lines below if you want to keep the longest protein isoform for each species
  
  unique_BLAST_results <- data.frame()
  spec <- unique(annotated_BLAST_results$critters) # list of all critters represented in the set
  for (i in 1:length(spec)) {
    # create a daraframe with all the sequences for one species
    BLAST_subset <- annotated_BLAST_results[annotated_BLAST_results$critters==spec[i],]
    # check that the sequence starts with a methionine
    BLAST_subset <- subset(BLAST_subset, grepl("^M", BLAST_subset$sequence, ignore.case=TRUE))
    # check whether the record still exists afterwards
    if (!is.na(BLAST_subset$sequence[1])) {
      if (keep_longest == TRUE) {
      # determine which of the list items is the longest sequence
      longest <- max(nchar(BLAST_subset$sequence[1:length(BLAST_subset$sequence)]))
      # write out the longest sequences into separate dataframe 
      # in case there are multiple sequences with the same (longest) length
      subset_longest <- BLAST_subset[nchar(BLAST_subset$sequence)==longest,] 
      unique_BLAST_results <- rbind(unique_BLAST_results,  subset_longest[1,])
      } else {
        # take only the first one, thus making sure there will be only one sequence per species
        unique_BLAST_results <- rbind(unique_BLAST_results,  BLAST_subset[1,])
      }
    }
  }
  # Usually, the original BLAST query sequence should be on top of the unique dataframe because that will have the best match in BLAST. 
  # However, it may not be the longest protein isoform in the organism and in that case the length-filter will have removed it in favor of a longer variant.
  # If required for downstream analysis, the simplest way is to add the sequence back to the start of the final fasta output file with a text editor.
  
  
  # reset the row names
  row.names(unique_BLAST_results) <- NULL
  
  # get the taxonomic ID's and if desired the common names where available
  unique_BLAST_results$species_uids <- as.numeric(get_uid(unique_BLAST_results$critters))
  # un-comment the line below if you also want to look up the common names (not required)
  #unique_BLAST_results$common_names <- as.character(sci2comm(unique_BLAST_results$critters, db = "ncbi"))
  
  # and write this to a file as dataframe
  filename_unique_BLAST_results <- paste0(names(protein_query_list[query_counter]),"_unique_BLAST_results.csv")
  write.csv(unique_BLAST_results, file=filename_unique_BLAST_results)
  
  #####################################################################################################################
  # At this point it is a good idea to verify if the sequence collection is as desired.                               #
  # In particular, there may be many more, closely related species than intended,                                     #
  # as this might limit the evolutionary "scope" present in the subsequent analysis.                                  #
  # It is possible to edit the collection by inserting / deleting rows into the .csv file with spreadsheet-software   #
  # and save it with the identical name.                                                                              #
  # To load the edited file: Just un-comment the line below and run the final steps of the script again.              #
  # Attention: If you indeed run a list of query proteins, this needs to be done explicitly and separately            #
  # for each candidate. Depending on how many candidates you are processing, an individual copy of the script for     #
  # each query sequence may be a better way to achieve your goals, rather than looping through a list.                #
  #####################################################################################################################
  
  # unique_BLAST_results <- read.csv(file=filename_unique_BLAST_result)   
  
  # create a fasta-file, as before with sink to avoid quotes
  filename_unique_BLAST_fasta <- paste0(names(protein_query_list[query_counter]),"_unique_BLAST_results.fa")
  sink(file=filename_unique_BLAST_fasta, append=FALSE, type="output", split = FALSE)
  for (i in 1:length(unique_BLAST_results$critters)){
    cat(paste(">",noquote(unique_BLAST_results$species_uids[i]),"|", noquote(unique_BLAST_results$accession[i]), "\n", sep=""))
    cat(paste(noquote(unique_BLAST_results$sequence[i]), "\n", sep=""))
  }
  sink(file = NULL)
  
  
  # Draw a phylogenetic tree of the species range covered by the analysis.
  # This is to obtain an overview of the "evolutionary reach" embraced by 
  # the recovered sequence set; it is NOT the phylogeny based on the BLAST hits!
  
 
  classifications <- classification(unique_BLAST_results$critters, db = "ncbi")
  
  lowest_common_rank <- lowest_common(unique_BLAST_results$species_uids, db= "ncbi", na.rm = TRUE)
  
  tree_info <- class2tree(classifications, check = FALSE) # check = FALSE avoids an error that can occur if only few species are considered
  
  filename_plot_linear <- paste0(names(protein_query_list[query_counter]),"_phylogeny.pdf")
  pdf(file=filename_plot_linear, width=10, height=10)
  plot(tree_info, cex = (40/length(unique_BLAST_results$critters)))
  title(main = paste("Query sequence:", names(protein_query_list[query_counter])))
  mtext(paste("lowest common rank:", lowest_common_rank$name, " - this tree is based on species phylogeny"), side=3, cex=0.8)
  dev.off()
  
  phylo_tree_obj <- tree_info$phylo
  phylo_tree_plot <- ggtree(phylo_tree_obj, layout="circular") + 
    geom_tiplab2(aes(angle=angle), size=240/length(unique_BLAST_results$critters)) + 
    theme_tree() +
    ggtitle(label= paste("Query sequence:", names(protein_query_list[query_counter])), subtitle = paste("lowest common rank:", lowest_common_rank, " - this tree is based on species phylogeny")) +
    theme(plot.title = element_text(hjust = 0.5)) +  # centered position
    theme(plot.subtitle = element_text(hjust = 0.5)) # centered position
  filename_plot_circle <- paste0(names(protein_query_list[query_counter]),"_circular_phylogeny.pdf")
  ggsave(filename=filename_plot_circle, 
         plot = phylo_tree_plot,
         device = pdf,
         width = 24,
         height = 30,
         units = "cm",
         dpi = 300)
  
  system.time({ Sys.sleep(120) }) # wait two minutes in between query sequences in the list to respect usage limits avoid being cut off from the NCBI BLAST service
}



