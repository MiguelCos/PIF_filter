---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
# Filtering of MaxQuant evidence.txt output by PIF for TMT quantitation

## Using the script:  

1. Open in an RStudio project folder.
2. Add the `evidence.txt` file in the same folder as the script.
3. Add the `msms.txt` file in the same folder as the script.
3. Open the script `pif_filtering_evidence_file_for_MSstats_fromMSMS.R` in RStudio.
4. Modify the Top_N criteria you want (line 9 of the script, 2 is default)
5. Modify the PIF cut-off you want to set (0.65 is default)
6. Execute the rest.

## Output:

A `.txt` will be generated in the same folder with the name `evidence_top_X_PIF_peptides.txt` where X would correspond to your desired Top_N criteria.  

__Note__: You might need to change the name of this file to `evidence.txt` when using it as input of other scripts for preprocessing. 

A secondary `filtered_out_proteins_non_quantifiable_peptides.txt` file will also be generated if the whole script is executed. This file would have the same format as the `evidence.txt` file and would only include the list of proteins that were excluded because were identified with unquantifiable peptides. 
