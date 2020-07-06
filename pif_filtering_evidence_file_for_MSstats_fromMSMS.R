### Script for PIF filtering on evidence.txt file from MaxQuant ####
### Miguel Cosenza 


### PLEASE SET THE FILTERING CRITERIA IN THE NEXT LINES #### 

# How many TOP_N peptides by PIF you want to keep by protein? (Default is 2)

top_pif_pepts <- 2

# Please set the minimum cutoff for the PIF values (Default is 0.65)

pif_cutoff <- 0.75

#### EXECUTION OF SCRIPT AND GENERATION OF FILTERED EVIDENCE FILE ####

# EXECUTE EVERY LINE CONSECUTIVELY 

### Install required packages if necessary ####

packages <- c("dplyr", "here", "stringr", "tidyr")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
      install.packages(setdiff(packages, rownames(installed.packages())))  
}

### Load required packages ####

library(dplyr)
library(stringr)
library(tidyr)
library(here)

### Load the msms.txt file ####

msms <- read.delim(file = here::here("msms.txt"),
                   sep = "\t",
                   stringsAsFactors = FALSE)

### Load evidence file ####

evidence <- read.delim(file = here::here("evidence.txt"),
                       sep = "\t",
                       stringsAsFactors = FALSE)

## Create the 'false' evidence file from the msms.txt file -----

# Get variable names present exclusively in the evidence file ----

only_in_evid <- setdiff(names(evidence),
                        names(msms)) 

only_in_evid <- c(only_in_evid, "id")

# Get mapping between MS.MS.IDs and ids in evidence file ----
# from selection of variable present only in evidence file 

# Get columns that are exclusive for the evidence file 
exclusive_evid <- dplyr::select(evidence,
                                Proteins,only_in_evid)

# Get msms file IDs associated to each peptide from the evidence file 
# Some peptides in the evidence file are selected from the several features in the msms file
# Here were are getting that info 
mmsms_id_split <- stringr::str_split(exclusive_evid$MS.MS.IDs,
                                     pattern = ";")

# The the maximum number of features that were used for a peptide from the msms file 
max_n_splits <- sapply(mmsms_id_split, length) %>% max()

# Empty string vector to deposit the name of the new variables created after separation by N features
new_vars <- as.character()

# Create a string vector with the names of the new variables after separation by N features
for (i in 1:max_n_splits){
      new_vars[i] <- paste0("msmsID_",i)
}


# Split the MSMSids column of the evidence file
msms_ids_cols_separ <- tidyr::separate(exclusive_evid,
                                       col = MS.MS.IDs,
                                       sep = ";",
                                       remove = FALSE,
                                       fill = "right",
                                       into = new_vars)  

# Get MSMSids info into an unique column after split  
evid_slim_wmsid <- tidyr::pivot_longer(msms_ids_cols_separ,
                            cols = new_vars,
                            names_to = "msmsids",
                            values_drop_na = TRUE,
                            values_to = "MS.MS.ID") %>% 
                   dplyr::select(-msmsids) %>%
                   dplyr::mutate(MS.MS.ID = str_trim(MS.MS.ID)) %>%
                   dplyr::mutate(MS.MS.ID = as.integer(MS.MS.ID))

# Rename MSMSid variable from the msms file to match the one in the evidence file for merging 
msms_mod <- dplyr::rename(msms, MS.MS.ID = id) %>%
            dplyr::select(-c(Proteins))


## Create 'disguised' evidence file containing all-features information from the msms file -----

disguised_evidence <- dplyr::left_join(msms_mod, evid_slim_wmsid,
                                     by = "MS.MS.ID") %>%
                      dplyr::select(names(evidence))

write.table(x = disguised_evidence,
            file = here::here("GBM_filtered_files/msms_disguised_as_evidence.txt"),
            sep = "\t",
            row.names = FALSE)

### Grouping and filtering msms file by PIF and Experiment -----

evidene_fil <- group_by(disguised_evidence, 
                        Experiment,
                        Proteins,
                        Raw.file) %>% # note to future Miguel: (sometimes it is better to not filter also by Raw.file, when there's only 1 experimental setting)
      dplyr::top_n(top_pif_pepts, wt = PIF) %>%
      dplyr::filter(PIF >= pif_cutoff) %>%
      ungroup()

# Save .txt filtered evidence file ####

write.table(x = evidene_fil,
            file = here::here(paste0("GBM_filtered_files/evidence_top_",as.character(top_pif_pepts),"pept_PIFfilter_",
                                     as.character(pif_cutoff),".txt")),
            sep = "\t",
            row.names = FALSE)

### Run this section of the script if you want to know the amount of proteins filtered out and which were those ####

# How many proteins were filtered out?
length(setdiff(evidence$Proteins,
               evidene_fil$Proteins))


# Which proteins were those? 

filtered_out_evidence <- dplyr::filter(evidence,
                                       Proteins %in% setdiff(evidence$Proteins,
                                                             evidene_fil$Proteins))

write.table(x = filtered_out_evidence,
            file = paste0("filtered_out_proteins_non_quantifiable_peptides.txt"),
            sep = "\t")