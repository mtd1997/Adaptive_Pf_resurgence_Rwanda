# Rename FASTQ files according to the sample barcode mapping table.

library(readr)

mapping <- read_csv("sample_barcode_mapping.csv")

for(i in seq_len(nrow(mapping))){

    old_name <- mapping$barcode[i]
    new_name <- paste0(mapping$SampleID[i], ".fastq.gz")

    if(file.exists(old_name)){

        file.rename(old_name, new_name)

        message(old_name, " -> ", new_name)

    }

}

message("Renaming completed.")
