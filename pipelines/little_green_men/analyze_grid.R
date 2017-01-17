#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(cowplot)

options("readr.num_columns" = 0)

top <- 7
total <- 15
p <- paste0("rs", 1:top)
n <- paste0("rs", (top+1):total)

parameters <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)


for (net in c("gs","gm","gi")){
  selectedSnps <- list()
  for (l.exp in parameters){
    for (e.exp in parameters){
      l <- paste0(1,paste0(rep("0",l.exp), collapse = ""))
      e <- paste0(1,paste0(rep("0",e.exp), collapse = ""))
      # SELECTED SNPS
      selectedSnpsFile <- paste0("gridsearch/BRCA_", net, ".l", l,".e", e,".scones.out.txt")
      this.selectedSnps <- read_tsv(selectedSnpsFile, comment = "#", col_names = F)
      if (nrow(this.selectedSnps)==0){

        tp <- NA
        tn <- NA
        
      } else {
        this.selectedSnps <- this.selectedSnps %>%
          set_colnames(c("snp","chr","pos"))
        
        tp <- length(intersect(this.selectedSnps$snp, p))
        fp <- length(intersect(this.selectedSnps$snp, n))
        tn <- length(n) - fp
        fn <- length(p) - tp
      }
      
      # TERMS
      termsFile <- paste0("gridsearch/BRCA_", net, ".l", l,".e", e,".scones.terms.txt")
      this.terms <- read_tsv(termsFile, col_names = F) %>%
        t %>% set_colnames(c("association","connectivity","sparsity"))

      selectedSnps[[paste0(l,e)]] <- data.frame(tpr = tp/length(p), tnr = tn/length(n),
                                                lambda = as.numeric(l), 
                                                eta = as.numeric(e),
                                                Association = this.terms[1],
                                                Connectivity = this.terms[2],
                                                Sparsity = this.terms[3])
      
    }
  }
  
  selectedSnps <- do.call("rbind",selectedSnps)
  
  tpr <- ggplot(selectedSnps, aes(lambda, eta, fill = tpr )) + 
    geom_tile(color = "gray") + 
    labs(x = "λ", y = "η", fill = "Sensitivity") +
    theme_minimal() +
    scale_x_log10() + 
    scale_y_log10()

  tnr <- ggplot(selectedSnps, aes(lambda, eta, fill = tnr )) + 
    geom_tile(color = "gray") + 
    labs(x = "λ", y = "η", fill = "Specificity") +
    theme_minimal() +
    scale_x_log10() + 
    scale_y_log10()
  
  terms <- selectedSnps %>%
    select(lambda, eta, Association, Connectivity, Sparsity) %>%
    gather(term, value, -lambda, -eta) %>%
    mutate(lambda = log10(lambda),
           term = recode(term, Association = "A", Connectivity = "C", Sparsity = "S")) %>%
    ggplot(aes(term, eta, fill = log10(value) )) +
    geom_tile(color = "gray") +
    labs(x = "Term", y = "η", fill = "log10(Term)", title = "log10(λ)") +
    theme_minimal() +
    scale_y_log10() +
    facet_grid(.~lambda) +
    theme(plot.title = element_text(hjust = 0.5))
  
  cv <- selectedSnps %>%
    group_by(lambda, eta) %>%
    summarise(cv = sd(c(Association, Connectivity, Sparsity))/mean(c(Association, Connectivity, Sparsity))) %>%
    ggplot(aes(lambda, eta, fill = cv )) + 
    geom_tile(color = "gray") + 
    labs(x = "λ", y = "η", fill = "CV") +
    theme_minimal() + 
    scale_y_log10() +
    scale_x_log10()
  
  plot_grid(tpr, tnr, terms, cv, 
            labels=c("Sensitivity","Specificity","Term order", "CV"), 
            label_size = 8) %>%
    ggsave(paste0(net, ".png"), ., width = 12, height = 9)
  
}