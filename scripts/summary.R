#!/usr/bin/env Rscript

library(magrittr)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

options("readr.num_columns" = 0)

args <- commandArgs(trailingOnly = TRUE)
# QUITAR
simTool <- "gametes"
software <- "beam"

quality <- list()

for (h in c("005","01","025","05","1","2")){
  for (maf in c("05","2","4")){
    for (N in c("20","100","1000")){
      # QUITAR
      for (modelNo in formatC(1:2, width = 2, flag = "0")){
        tests <- data_frame()
        # QUITAR
        for (repNo in formatC(1:3, width = 3, flag = "0")){
          # PLINK
          ###########################
          if (software == "plink"){
            replicate <- paste0("h",h,"_maf",maf,"_N",N,"_EDM-",modelNo,"_",repNo,".plink.txt.epi.cc")
            replicatePath <- paste0("populations/",simTool,"/plink/",replicate)

            if (file.exists(replicatePath)){
              # read the tests, get the pvalues
              thisSample <- read_fwf(replicatePath, fwf_widths(c(4,5,5,5,13,13,13)), skip = 1) %>%
                set_colnames(c("chr1", "snp1", "chr2", "snp2", "beta", "stat", "p")) %>%
                select(snp1,snp2,p) %>%
                mutate(., p = p * nrow(.))
            }
          }
          # BEAM
          ###########################
          else if (software == "beam"){
            replicate <- paste0("h",h,"_maf",maf,"_N",N,"_EDM-",modelNo,"_",repNo,".beam.txt")
            replicatePath <- paste0("populations/",simTool,"/beam/",replicate)

            if (file.exists(replicatePath) & file.info(replicatePath)$size > 0){
              # read the tests, get the pvalues
              thisSample <- read_tsv(replicatePath, col_types = cols(Marker_Group = "c", Pvalue = "d")) %>%
                mutate(Marker_Group = gsub(",$", "", Marker_Group)) %>%
                separate(Marker_Group, c("snp1", "snp2", "snp3"), ",") %>%
                filter(is.na(snp3)) %>%
                select(snp1, snp2, Pvalue) %>%
                rename(p = Pvalue) %>%
                mutate(., p = p * nrow(.))
            }
          }

          # create table
          if (nrow(tests) == 0){
            tests <- thisSample
          } else {
            tests <- merge(tests,thisSample,by = c("snp1","snp2"), suffixes = c("",repNo))
          }
        }

        if (nrow(tests) > 0){
          pvals <- tests %>% select(starts_with("p")) %>% as.matrix
          models <- tests %>%
            select(-starts_with("p")) %>%
            mutate(NumPositives = rowSums(pvals < 0.05, na.rm = T))

          # quality measures
          if (software == "plink"){
            causalSnpRow <- models$snp1 == "M0P0" & models$snp2 == "M0P1"
          } else if (software == "beam"){
            M0P0 <- as.character(as.numeric(N) - 2)
            M0P1 <- as.character(as.numeric(N) - 1)
            causalSnpRow <- models$snp1 %in% c(M0P0,M0P1) & models$snp2 %in% c(M0P0,M0P1)
          }

          TP <- ifelse(sum(causalSnpRow), models$NumPositives[causalSnpRow], 0)
          FN <- ncol(pvals) - TP
          FP <- sum(models$NumPositives[! causalSnpRow])
          TN <- (nrow(pvals) - 1) * ncol(pvals) - FP

          acc <- (TP + TN)/(TP + TN + FP + FN)
          tpr <- TP/(TP + FN)
          tnr <- TN/(TN + FP)

          quality[[paste(h, maf, modelNo, N)]] <- data.frame( h = h, maf = maf, model = modelNo, N = N,
                                                              accuracy = acc, sensitivity = tpr, specificity = tnr)
        }
      }
    }
  }
}

quality.sum <- do.call("rbind", quality) %>%
  as.data.frame %>%
  group_by(h, maf, model, N) %>%
  summarise(acc = median(accuracy), min_acc = min(accuracy), max_acc = max(accuracy),
            tpr = median(sensitivity), min_tpr = min(sensitivity), max_tpr = max(sensitivity),
            tnr = median(specificity), min_tnr = min(specificity), max_tnr = max(specificity)) %>%
  ungroup %>%
  mutate(h = paste0("h = 0.", h), maf = paste0("MAF = 0.", maf))

ggplot(quality.sum, aes(x = N, y = tpr, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_tpr, ymin=min_tpr), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Sensitivity", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/", paste(software,simTool,"sensitivity.png", sep=".")))

ggplot(quality.sum, aes(x = N, y = tnr, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_tnr, ymin=min_tnr), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Specificity", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/", paste(software,simTool,"specificity.png", sep=".")))

ggplot(quality.sum, aes(x = N, y = acc, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_acc, ymin=min_acc), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Accuracy", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/", paste(software,simTool,"accuracy.png", sep=".")))


# quality.sum <- do.call("rbind", quality) %>%
#   as.data.frame %>%
#   group_by(h, maf, N) %>%
#   summarise(acc = median(accuracy), sd_acc = sd(accuracy),
#             tpr = median(sensitivity), sd_tpr = sd(sensitivity),
#             tnr = median(specificity), sd_tnr = sd(specificity)) %>%
#   ungroup %>%
#   mutate(h = paste0("h = 0.", h), maf = paste0("MAF = 0.", maf))
