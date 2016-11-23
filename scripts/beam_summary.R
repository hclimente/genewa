#!/usr/bin/env Rscript

library(magrittr)
library(readr)
library(dplyr)
library(ggplot2)

options("readr.num_columns" = 0)

args <- commandArgs(trailingOnly = TRUE)
simTool <- args[1]

beam.quality <- list()

for (h in c("005","01","025","05","1","2")){
  for (maf in c("05","2","4")){
    for (N in c("20","100","1000")){
      for (modelNo in formatC(1:10, width = 2, flag = "0")){
        beam.tests <- data_frame()

        for (repNo in formatC(1:100, width = 3, flag = "0")){
          replicate <- paste0("h",h,"_maf",maf,"_N",N,"_EDM-",modelNo,"_",repNo,".beam.txt")
          replicatePath <- paste0("populations/",simTool,"/beam/",replicate)

          if (file.exists(replicatePath)){
            # read the tests, get the pvalues
            thisSample <- read_tsv(replicatePath) %>%
              mutate(Marker_Group = gsub(",$", "", Marker_Group)) %>%
              separate(Marker_Group, c("snp1", "snp2", "snp3"), ",") %>%
              filter(is.na(snp3)) %>%
              select(snp1, snp2, Pvalue)

            # create table
            if (nrow(beam.tests) == 0){
              beam.tests <- thisSample
            } else {
              beam.tests <- merge(beam.tests,thisSample,by = c("snp1","snp2"), suffixes = c("",repNo))
            }
          }
        }

        if (nrow(beam.tests) > 0){
          pvals <- beam.tests %>% select(starts_with("Pvalue")) %>% as.matrix
          models <- beam.tests %>%
            select(-starts_with("Pvalue")) %>%
            mutate(NumPositives = rowSums(pvals < 0.05, na.rm = T))

          # count cases with NAs (no test performed)
          M0P0 <- as.character(N - 2)
          M0P1 <- as.character(N - 1)
          naPositives <- sum(is.na(pvals[models$snp1 == M0P0 & models$snp2 == M0P1]), na.rm = T)
          naNegatives <- pvals[models$snp1 != M0P0 | models$snp2 != M0P1,] %>%
            is.na %>%
            rowSums %>%
            sum

          # quality measures
          TP <- models$NumPositives[models$snp1 == M0P0 & models$snp2 == M0P1] - naPositives
          FN <- ncol(pvals) - models$NumPositives[models$snp1 == M0P0 & models$snp2 == M0P1] - naPositives
          FP <- sum(models$NumPositives[models$snp1 != M0P0 | models$snp2 != M0P1]) - naNegatives
          TN <- (nrow(pvals) - 1) * ncol(pvals) - naNegatives

          acc <- (TP + TN)/(TP + TN + FP + FN)
          tpr <- TP/(TP + FN)
          tnr <- TN/(TN + FP)

          beam.quality[[paste(h, maf, modelNo, N)]] <- data.frame( h = h, maf = maf, model = modelNo, N = N,
                                                                  accuracy = acc, sensitivity = tpr, specificity = tnr)
        }
      }
    }
  }
}

beam.quality.sum <- do.call("rbind", beam.quality) %>%
  as.data.frame %>%
  group_by(h, maf, model, N) %>%
  summarise(acc = median(accuracy), min_acc = min(accuracy), max_acc = max(accuracy),
            tpr = median(sensitivity), min_tpr = min(sensitivity), max_tpr = max(sensitivity),
            tnr = median(specificity), min_tnr = min(specificity), max_tnr = max(specificity)) %>%
  ungroup %>%
  mutate(h = paste0("h = 0.", h), maf = paste0("MAF = 0.", maf))

ggplot(beam.quality.sum, aes(x = N, y = tpr, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_tpr, ymin=min_tpr), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Sensitivity", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/beam.",simTool,".sensitivity.png"))

ggplot(beam.quality.sum, aes(x = N, y = tnr, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_tnr, ymin=min_tnr), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Specificity", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/beam.",simTool,".specificity.png"))

ggplot(beam.quality.sum, aes(x = N, y = acc, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = max_acc, ymin=min_acc), position = "dodge") +
  ylim(c(0,1)) +
  labs(x = "# SNPs", y = "Accuracy", fill = "Model") +
  facet_grid(h ~ maf) +
  theme_minimal()
ggsave(paste0("results/sota_benchmark/beam.",simTool,".accuracy.png"))
