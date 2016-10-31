source("~/libs/wisdom/r/data_analysis_environment.R")

args <- commandArgs(trailingOnly = TRUE)
ped <- args[1]
beam <- args[2]

snps <- read_tsv(ped, col_names = F) %>% 
  # remove columns related to family, individual, father, mother and sex
  select(-X1, -X2, -X3, -X4, -X5)

# fetch phenotype information and remove that column
phenotypes <- snps$X6 %>% rep(each=2)
snps <- snps %>% select(-X6)

# take separately the two sets of chromosomes
haplo1 <- snps[,c(T,F)] %>% t
colnames(haplo1) <- paste0("h1",1:dim(haplo1)[2])
haplo2 <- snps[,c(F,T)] %>% t
colnames(haplo2) <- paste0("h2",1:dim(haplo2)[2])

# combine the data to create beam format
intertwinedCols <- paste0(c("h1","h2"),rep(1:dim(haplo2)[2], each=2))
cbind(haplo1, haplo2) %>%
  .[,intertwinedCols] %>%
  rbind(phenotypes, .) %>%
  write_tsv(beam)