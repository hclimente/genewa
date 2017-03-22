#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(cowplot)

options("readr.num_columns" = 0)

args <- commandArgs(trailingOnly = TRUE)
top <- as.numeric(args[1])
total <- as.numeric(args[2])
net <- args[3]
lambda_base <- as.numeric(args[4])
eta_base <- as.numeric(args[5])

p <- paste0("rs", 1:top)
n <- paste0("rs", (top+1):total)

parameters <- seq(-5,5)

##############################
#  Read individual analyses  #
##############################

selectedSnps <- lapply(parameters, function(l){
  lapply(parameters, function(e){
    # Selected SNPs
    this.selectedSnps <- paste0("BRCA_", net, ".l", lambda_base,"e", l,".e", eta_base, "e", e,".scones.out.ext.txt") %>%
      read_tsv(comment = "#", col_names = F, col_types = cols(X5 = col_double()))
    if (nrow(this.selectedSnps)==0){
      
      tp <- NA
      tn <- NA
      
    } else {
      this.selectedSnps <- this.selectedSnps %>%
        set_colnames(c("snp","chr","pos","selected","skat")) %>%
        filter(selected >= 1)
      
      tp <- length(intersect(this.selectedSnps$snp, p))
      fp <- length(intersect(this.selectedSnps$snp, n))
      tn <- length(n) - fp
      fn <- length(p) - tp
      
    }
    
    # Terms
    this.terms <- paste0("BRCA_", net, ".l", lambda_base,"e", l,".e", eta_base, "e", e,".scones.terms.txt") %>%
      read_tsv(col_names = F, col_types = cols(X1 = col_double())) %>%
      t %>% set_colnames(c("association","connectivity","sparsity"))
    
    data.frame(tpr = tp/length(p), 
               tnr = tn/length(n),
               lambda = lambda_base**(l),
               eta = eta_base**(e),
               Association = this.terms[1],
               Connectivity = this.terms[2],
               Sparsity = this.terms[3],
               proportionPositives = (tp + fp)/total)

    
  }) %>% do.call("rbind",.)
}) %>% do.call("rbind",.)

l <- log(selectedSnps$lambda[which.max((selectedSnps$tpr + selectedSnps$tnr)/2)], base = lambda_base)
e <- log(selectedSnps$eta[which.max((selectedSnps$tpr + selectedSnps$tnr)/2)], base = eta_base)

bestSnps <- paste0("BRCA_", net, ".l", lambda_base,"e", l,".e", eta_base, "e", e,".scones.out.ext.txt") %>%
  read_tsv(comment = "#", col_names = F, col_types = cols(X5 = col_double())) %>%
  set_colnames(c("snp","chr","pos","selected","skat")) %>%
  filter(selected >= 1) %>%
  .$snp

###########################
#    Read SKAT scores     #
###########################
scores <- paste0("BRCA_", net, ".l", lambda_base,"e1.e", eta_base, "e1.scones.out.ext.txt") %>%
  read_tsv(comment = "#", col_names = F, col_types = cols(X5 = col_double()))  %>%
  set_colnames(c("snp","chr","pos","selected","skat")) %>%
  mutate(causal = ifelse(snp %in% p, "Yes", "No"))

###########################
#    Read interactions    #
###########################
L <- paste0("BRCA_", net, ".scones.L.txt") %>%
  read_tsv %>%
  apply(1, function(x) x > 0) %>%
  as.data.frame %>%
  set_colnames(., sub("V","rs", colnames(.))) %>%
  mutate(., rs = colnames(.))

L %>%
  gather(rs) %>%
  set_colnames(c("rs1", "rs2", "interaction")) %>%
  merge(scores %>% select(snp,chr,pos), by.x="rs1", by.y="snp") %>%
  merge(scores %>% select(snp,chr,pos), by.x="rs2", by.y="snp", suffixes = c("1", "2")) %>%
  arrange(chr1,pos1,chr2,pos2) %>%
  mutate(rs1 = factor(rs1, levels = scores %>% arrange(chr,pos) %>% .$snp),
         rs2 = factor(rs2, levels = scores %>% arrange(chr,pos) %>% .$snp),
         interaction = ifelse(interaction, "Yes", "No")) %>%
  ggplot(aes(x=rs1, y=rs2, fill=as.character(chr1), alpha=interaction)) +
    geom_tile() +
    # scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank()) + 
    labs(fill = "Chromosome", alpha = "Interaction")
ggsave(paste0(net,".interactions.png"))

###########################
#          Plots          #
###########################
tpr <- ggplot(selectedSnps, aes(lambda, eta, fill = sqrt(tpr) )) +
  geom_tile(color = "gray") +
  labs(x = "λ", y = "η", fill = "sqrt(sensitivity)") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()

tnr <- ggplot(selectedSnps, aes(lambda, eta, fill = sqrt(tnr) )) +
  geom_tile(color = "gray") +
  labs(x = "λ", y = "η", fill = "sqrt(specificity)") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()

terms <- selectedSnps %>%
  select(lambda, eta, Association, Connectivity, Sparsity) %>%
  gather(term, value, -lambda, -eta) %>%
  mutate(lambda = log(lambda, base = lambda_base),
         term = recode(term, Association = "A", Connectivity = "C", Sparsity = "S")) %>%
  ggplot(aes(term, eta, fill = log10(value) )) +
  geom_tile(color = "gray") +
  labs(x = "Term", y = "η", fill = "log10(Term)", title = paste0("log", lambda_base,"(λ)")) +
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

numSnps <- ggplot(selectedSnps, aes(lambda, eta, fill = sqrt(proportionPositives) )) +
  geom_tile(color = "gray") +
  labs(x = "λ", y = "η", fill = "Positive SNPs (sqrt(%))") +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10()

skat <- scores %>%
  arrange(chr, pos) %>%
  mutate(gPos = paste(chr, pos),
         gPos = factor(gPos, levels = gPos),
         selected = ifelse(snp %in% bestSnps, "Yes", "No")) %>%
  ggplot(aes(gPos, skat, color = as.character(chr), alpha = causal, shape = selected )) +
  geom_point() +
  labs(x = "Genomic position", y = "Score", alpha = "Causal", 
       color = "Chromosome", shape = "Selected") +
  theme_minimal() + 
  theme(axis.text.x=element_blank()) 

plot_grid(tpr, tnr, numSnps, terms, cv, skat,
          labels=c("Sensitivity","Specificity", "Positives", "Term order", "CV", 
                   paste0("Score at η = ", eta_base,"e",e," λ = ", lambda_base,"e",l)),
          label_size = 8, nrow = 2) %>%
  ggsave(paste(net,"png", sep = "."), ., width = 18, height = 9)