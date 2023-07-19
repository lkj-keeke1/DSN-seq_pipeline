library(tidyverse)

# input

## merge the sites info of filtered.table
edit_sites <- c()
for (sample in c("Ath_DSN_merge.")) {
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

for (sample in c("Ath_DSN_merge.")) {
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

save(edit_all, file = "Ath_edit_all.Rdata")
