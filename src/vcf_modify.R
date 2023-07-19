args <- commandArgs(trailingOnly = T)


vcf <- args[1]                      # Varscan2 output vcf file
reference_length_table <- args[2]   # gene length table in the multi-FASTA genome, including flanking sequences
outtable_prefix <- args[3]          # output prefix

f_p_value <- as.numeric(args[4])   # filtering standard of p-value (calculated by Varscan2)
f_freq <- as.numeric(args[5])     # filtering standard of frequency/editing extent
f_depth <- as.numeric(args[6])     # filtering standard of sequencing depth


library(dplyr)
library(tidyr)


########## Step1: Import Data ##########

# import vcf table
vcftable  <- read.table(vcf, comment.char = "#", sep = "\t", header = F) %>% 
  separate(col = V10, into = c("GT", "GQ", "SDP", "DP", "RD", "AD", "FREQ", "PVAL", "RBQ", "ABQ", "RDF", "RDR", "ADF", "ADR"),sep = ":") %>% 
  separate(col = FREQ, into = c("freq", NA), sep = "%") %>% 
  mutate(gene = V1,
         posi = as.numeric(V2) - 100, # 相对 ATG 的位置
         depth = as.numeric(DP),
         ref_num = as.numeric(RD),
         var_num = as.numeric(AD),
         edit = paste(V4, V5, sep = ""),
         freq = as.numeric(freq),
         pvalue = as.numeric(PVAL),
         edit_id = paste(gene, posi, edit, sep = " ")
         ) %>% 
  select(edit_id,
         gene,
         posi,
         edit,
         depth,
         ref_num,
         var_num ,
         freq,
         pvalue)

# import length table 
length_table <- read.table(reference_length_table,
                           header = F,
                           sep = "\t",
                           col.names = c("gene", "length"))


########## Step2: Filtering ##########

# calculate the length of CDSs 
vcftable <- left_join(vcftable, length_table, by = "gene") %>% 
  mutate(length = length - 200) # `CDS length` = `full length in multi-FASTA file` - `flanking sequences length`

# filtering step
vcftable_filter <- filter(vcftable, 
                          edit == "CT",            # only C-to-T editing sites
                          posi >0, posi <= length, # editing sites on CDSs
                          pvalue <= f_p_value,     # p-value
                          freq >= f_freq,          # frequency/editing extent
                          depth >= f_depth         # sequencing depth
                          )

# unqualified sites
vcf_filtered_out <- filter(vcftable, !(edit_id %in% vcftable_filter$edit_id))


########## Step2: Output ##########

filtered_table <- paste(outtable_prefix, ".filtered.table", sep = "")  # qualified sites 
filtered_out_table <- paste(outtable_prefix, ".out.table", sep = "")   # unqualififed sites

write.table(file = filtered_table,
            x = vcftable_filter,
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)

write.table(file = filtered_out_table,
            x = vcf_filtered_out,
            sep = "\t",
            col.names = T,
            row.names = F,
            quote = F)
