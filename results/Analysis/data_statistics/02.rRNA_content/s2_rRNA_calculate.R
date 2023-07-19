rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)


colnames <- c('sample',
              'total_reads',
             'secondary',
             'suplementary',
             'duplicates',
             'mapped_reads',
             'paired_reads',
             'read1',
             'read2',
             'properly_paired_reads',
             'R1+R2_mapped',
             'singletons',
             'different_chr_mapping',
             'different_chr_mapping(mapQ>=5)')


# input merged flagstat table
rRNA_c <- as.data.frame(t(read.table("rRNA_content.summary",
                                           header = F)))
colnames(rRNA_c) <- colnames

# calculate rRNA content
rRNA_c <- mutate(rRNA_c, 
                 mapped_reads = as.numeric(mapped_reads),
                 total_reads = as.numeric(total_reads),
                 rRNA_content = (mapped_reads / total_reads) * 100) %>% 
  select(sample, rRNA_content)

# output table
write.table(x = rRNA_c,
            file = "rRNA_content.table",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)