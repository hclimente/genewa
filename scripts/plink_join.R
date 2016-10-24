library(magrittr)
library(readr)
library(dplyr)

plink <- list()

for (h in c("005","01","025","05","1","2")){
  for (maf in c("2","4")){
    for (modelNo in formatC(1:10, width = 2, flag = "0")){
      for (N in c("20","100","1000")){
        model <- paste0("h",h,"_maf",maf,"_N",N,"_EDM-",modelNo)
        for (repNo in formatC(1:100, width = 3, flag = "0")){
          replicate <- paste0(model,"_",repNo,".txt.ped.plink.epi.qt")
          
          plink[[replicate]] <- read_fwf(paste(model,replicate, sep = "/"), fwf_widths(c(4,5,5,5,13,13,13)), skip = 1) %>% 
            set_colnames(c("chr1", "snp1", "chr2", "snp2", "beta", "stat", "p")) %>% 
            mutate(h = h, maf = maf, N = 100, model = modelNo, sample = repNo)
        }
      }
    }
  }
}

plink <- do.call("rbind", plink)

write_tsv(plink, "alltogether.tsv")