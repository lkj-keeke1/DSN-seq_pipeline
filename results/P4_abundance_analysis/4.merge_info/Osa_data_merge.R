

library(tidyverse)

# input

## new classical sites

classical_sites_file <- "~/projects/DSN-seq/data/ref/Osa/Osa_new_canonical_sites.txt"

classcial_sites <- read.table(classical_sites_file, header = T, sep = "\t") %>% 
  mutate(chr = if_else(Organelle == "Mt", "_mt", "_pt"),
         edit_id = paste(Gene, chr, sep = ""),
         edit_id = paste(edit_id, `Position.in.CDS`, "CT", sep = " ")) %>% pull(edit_id)


## merge the sites info of filtered.table
edit_sites <- c()
for (sample in c("Osa_DSN_merge.", "Osa_Ribo-off_merge.")) {
  for (seed in c("11", "22", "33")) {
    for (frac in 1:9) {
      input_file <- paste("../3.filter/", sample, seed, "_", frac, ".filtered.table", sep = "")
      tmp_table <- read.table(input_file, sep = "\t", header = T) %>% pull(edit_id)
      edit_sites <- c(edit_sites,tmp_table)
    }
  }
}

edit_all <- distinct(tibble(edit_sites)) %>% 
  dplyr::rename(edit_id = edit_sites)


# data process

## annotate classical sites

edit_all <- mutate(edit_all, classical = edit_all$edit_id %in% classcial_sites)


for (sample in c("Osa_DSN_merge.", "Osa_Ribo-off_merge.")) {
  for (seed in c("11", "22", "33")) {
    for (frac in 1:9) {
      input_file <- paste("../3.filter/", sample, seed, "_", frac, ".filtered.table", sep = "")
      tmp_id <- read.table(input_file, sep = "\t", header = T) %>%pull(edit_id)
      edit_tmp <- select(edit_all,edit_id) %>% mutate(tmp = edit_id %in% tmp_id)
      colnames(edit_tmp) <- c("edit_id", paste(sample, seed, "_", frac, sep = ""))
      edit_all <- left_join(edit_all, edit_tmp, by = "edit_id")
      }
    }
  }



# output 

save(edit_all, file = "Osa_edit_all.Rdata")
