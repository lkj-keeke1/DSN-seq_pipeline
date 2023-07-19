##### annotate amino acid change #####

rm(list=ls())
options(stringsAsFactors = F)


library(tidyverse)
library(readr)
library(IRanges)

##------------------- function ------------------------

codon_edit <- function(pri_codon, edited_posi, orginal_base, final_base) {
  
  # change the pri_codon (primary) to aft_codon (after editing)
  
  # Parameter:
  ## pri_codon: str, primary codon sequence of edited codon
  ## edited_posi: numeric, the position of edited base on the codon
  ## original_base: str, theoretically orgnical base on the edited position of codon
  ## final_base: str, the final base after editing on the edited positon of codon
  
  # Return：
  ## edited codon
  
  if (substring(pri_codon, edited_posi, edited_posi) == orginal_base) {
    substring(pri_codon, edited_posi, edited_posi) <- final_base
  } else {
    # if the original base != the corresponding base on the codon: wrong genome sequence or wrong editing site
    print(paste("Not consistent original_base: \nref: ", orginal_base, "\nedited_posi: ", edited_posi, "\npri_codon: ", pri_codon, sep = ""))
    return("wrong_record")  
  }
  
  
  return(pri_codon)
  
}

##---------------- Input -----------------------

load("tmp1_gene_info_annotated.Rdata")

# genetic code table
## e.g.：TTT Phe F
genetic_code_file <- "/home/hj-z2/projects/DSN-seq/data/2022-07-29_genetic_code_tbl/genetic_code.txt"
genetic_code_table <- read.delim(genetic_code_file, header=FALSE)

# primary codon sequences of edited sites
codonseq_file <- "tmp_pre_codon.table"
editing_pri_codon <- read.table(codonseq_file, header = F, sep = "\t") %>% 
  select(V2,V3)
colnames(editing_pri_codon) <- c('edit_id', 'pri_codon')


##---------------- data process -----------------------

if (TRUE) {
  
  # merge the pri_codon 
  all_sites <- left_join(all_sites, editing_pri_codon, by = "edit_id")
  
  # initialize amino acid change-related column
  all_sites <- mutate(all_sites,
                      pri_amino_acid = NA, # amino acid translated from primary codon
                      aft_codon = NA,      # edited codon sequence
                      aft_amino_acid = NA  # amino acid translated from edited codon
  )
  
  # obtain pri_amino_acid
  for (i in 1:nrow(all_sites)) {
    if (!is.na(all_sites$pri_codon[i]) && all_sites$posi[i] > 0) {
      all_sites$pri_amino_acid[i] <- genetic_code_table$V3[which(genetic_code_table$V1 == all_sites$pri_codon[i])]
    }
  }
  
  # arrange the editing sites based on seq and position
  ## this will put adjacent sites together
  ## one codon can contain >1 editing sites, all sites will be used for annotating amino acid change
  all_sites <- arrange(all_sites, seq, posi)
  
  
  # obtain aft_codon and aft_amino_acid

  i= 1 # row number

  while (i <= nrow(all_sites)) {
    
    if (!is.na(all_sites$pri_amino_acid[i])) {
      
      # save row numbers of editing sites in one codon
      is <- c(i)
      # unique ID of each codon, used to determine adjacent editing sites in one codon 
      i_codon_id <- paste(all_sites$seq[i],all_sites$codon_num[i], sep = '')

      # read the codon information
      i_pri_codon <- all_sites$pri_codon[i]                     # primary codon sequence
      codon_edited_posi <- as.numeric(all_sites$codon_posi[i])  # edited position
      
      pri_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][1] # primary base
      aft_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][2] # base after editing
      
      # edit the codon sequence
      i_pri_codon <- codon_edit(i_pri_codon,codon_edited_posi, pri_base, aft_base)
      
      
      # determine adjacent editing sites in the same codon
      ## determine whether the i+1 row-recored editing site is still in the same codon
      next_i_codon_id <- paste(all_sites$seq[i+1],all_sites$codon_num[i+1], sep = '')
      if (next_i_codon_id == i_codon_id) {
        
        i <- i + 1
        is <- c(is, i)
        
        codon_edited_posi <- as.numeric(all_sites$codon_posi[i])
        
        pri_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][1]
        aft_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][2]
        
        # edit the codon sequence again
        i_pri_codon <- codon_edit(i_pri_codon, codon_edited_posi, pri_base, aft_base)
        
        
        ## determine whether the i+2 row-recored editing site is still in the same codon: this is corresponding to CCC -> TTT 
        next_i_codon_id <- paste(all_sites$seq[i+1],all_sites$codon_num[i+1], sep = '')
        if (next_i_codon_id == i_codon_id) {
          i <- i + 1
          is <- c(is, i)
          
          codon_edited_posi <- as.numeric(all_sites$codon_posi[i])
          
          pri_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][1]
          aft_base = str_split(all_sites$editing_type[i], pattern = "")[[1]][2]
          
          # edit the codon sequence for the third time
          i_pri_codon <- codon_edit(i_pri_codon,codon_edited_posi, pri_base, aft_base)
        }
      }
      
      
      # annotate the edited codon for all the editing sites in the same codon(including only 1 sites in one codon)
      for (j in is) {
        all_sites$aft_codon[j] <- i_pri_codon
      }
      
      
      i = i + 1
      
    } else {
      
      # this editing site not belong to protein-coding genes
      i = i + 1
      
    }
    
  }
  
  # obtain the aft_amino_acid based on aft_codon
  ## in the meantime, annotate the codon_change and aa_change column
  all_sites <- mutate(all_sites,
                      codon_change = NA, # codon change: pri_codon -> aft_codon
                      aa_change = NA,    # amino acid change: pri_amino_acid -> aft_amino_acid
                      syno = NA)         # whether synonymous amino acid change (amino acid does not change)
  
  for (i in 1:nrow(all_sites)) {
    if (!is.na(all_sites$aft_codon[i]) &&
        all_sites$aft_codon[i] != "wrong_record") {
      
      pri_codon <- all_sites$pri_codon[i]
      aft_codon <- all_sites$aft_codon[i]
      
      pri_aa <- all_sites$pri_amino_acid[i]
      aft_aa <- genetic_code_table$V3[which(genetic_code_table$V1 == all_sites$aft_codon[i])]
      
      all_sites$aft_amino_acid[i] <- aft_aa
      
      all_sites$codon_change[i] <- paste(pri_codon, "->", aft_codon, sep = "")
      all_sites$aa_change[i] <- paste(pri_aa, "->", aft_aa, sep = "")
      
      if (pri_aa == aft_aa) {
        all_sites$syno[i] <- "syn"
        
      } else {
        all_sites$syno[i] <- "non-syn"
      }
    }
  }
  
}




# rename
all_sites_annotation <- all_sites

##---------------- output -----------------------

save(all_sites_annotation,
     file = "tmp2_aa_change_annotated.Rdata")


