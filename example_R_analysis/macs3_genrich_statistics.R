rm(list = ls())
library(edgeR)
library(readr)
library(stringr)
library(pheatmap)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ggplot2)
library(clusterProfiler)
library(dplyr)



load("./analysis/intermediate_dataobjects.Rdat")
###### model fit;

##### First set of statistics does the celltype effect

design <- model.matrix(~ celltype, data = dge_macs3$samples)

macs3_fit = glmQLFit(dge_macs3, design = design)
macs3_res = glmQLFTest(macs3_fit, coef = 2)
genrich_fit = glmQLFit(dge_genrich, design = design)
genrich_res = glmQLFTest(genrich_fit, coef = 2)

macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0) %>%
  pull(logFC_dir) %>%
  table() # 199 False, 147 TRue

genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0) %>%
  pull(logFC_dir) %>%
  table() # 137 false, 188 true



data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, "")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)
ggsave("./plots_and_tables/volcano_macs3.png", dpi = 300,
       width = 6, height = 4)

library(tidyr)

data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(annotation = ifelse(FDR < 0.05 & abs(logFC) > 1,
                             annotation, NA)) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL,
             color = annotation))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)
ggsave("./plots_and_tables/volcano_annotations_macs3.png", dpi = 300,
       width = 7, height = 4)



data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, "")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)
ggsave("./plots_and_tables/volcano_genrich.png", dpi = 300,
       width = 6, height = 4)


data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(annotation = ifelse(FDR < 0.05 & abs(logFC) > 1,
                             annotation, NA)) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL,
             color = annotation))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)
ggsave("./plots_and_tables/volcano_annotations_genrich.png", dpi = 300,
       width = 7, height = 4)

### bar charts of peak location annotations
data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(abs(logFC) > 1 & FDR < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "higher in ST", "higher in SR")) %>%
  count(annotation, direction) %>%
  mutate(n = ifelse(direction == "higher in SR", -n, n)) %>%
  ggplot(aes(x = annotation, y = n, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("higher in SR" = "#D55E00", "higher in ST" = "#0072B2")) +
  theme_bw(10) +
  labs(x = "Annotation", y = "# significant peaks", fill = "binding") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
ggsave("./plots_and_tables/barchart_annotations_macs3.png", dpi = 300,
       width = 5, height = 2.5)



data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(abs(logFC) > 1 & FDR < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "higher in ST", "higher in SR")) %>%
  count(annotation, direction) %>%
  mutate(n = ifelse(direction == "higher in SR", -n, n)) %>%
  ggplot(aes(x = annotation, y = n, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("higher in SR" = "#D55E00", "higher in ST" = "#0072B2")) +
  theme_bw(10) +
  labs(x = "Annotation", y = "# significant peaks", fill = "binding") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
ggsave("./plots_and_tables/barchart_annotations_genrich.png", dpi = 300,
       width = 5, height = 2.5)


###### some point plots for top genes for gut check
cpm_macs3 = cpm(dge_macs3, log = TRUE)
cpm_genrich = cpm(dge_genrich, log = TRUE)
macs3_sig = macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0)  %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)

cpm_macs3[c(macs3_sig$peak[macs3_sig$logFC > 0][1:20],
            macs3_sig$peak[macs3_sig$logFC < 0][1:20]),] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "CPM") %>%
  left_join(metadata_sub) %>%
  left_join(macs3_sig) %>%
  separate(annotation, into = c("annotation"), sep = "\\(", extra = "drop") %>%
  arrange(logFC) %>%
  mutate(annotation = factor(annotation, levels = unique(annotation)),
         SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)),
         peak = factor(peak, levels = unique(peak))) %>%
  ggplot(aes(x = celltype, color = drug, y = CPM))+
  geom_point(size = 1.5,position = position_dodge(width = 0.2))+
  facet_wrap(peak ~ SYMBOL ~ annotation, scales = "free_y", nrow =4)+
  theme_bw(5.5)+
  labs(y = "log10(CPM)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(3, "mm"))
ggsave("./plots_and_tables/cpm_plots_top_DEGs_macs3.png", dpi = 300,
       width = 8, height = 4)

genrich_sig = genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0)  %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)

cpm_genrich[c(genrich_sig$peak[genrich_sig$logFC > 0][1:20],
              genrich_sig$peak[genrich_sig$logFC < 0][1:20]),] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "CPM") %>%
  left_join(metadata_sub) %>%
  left_join(genrich_sig) %>%
  separate(annotation, into = c("annotation"), sep = "\\(", extra = "drop") %>%
  arrange(logFC) %>%
  mutate(annotation = factor(annotation, levels = unique(annotation)),
         SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)),
         peak = factor(peak, levels = unique(peak))) %>%
  ggplot(aes(x = celltype, color = drug, y = CPM))+
  geom_point(size = 1.5,position = position_dodge(width = 0.2))+
  facet_wrap(peak ~ SYMBOL ~ annotation, scales = "free_y", nrow =4)+
  theme_bw(5.5)+
  labs(y = "log10(CPM)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(3, "mm"))
ggsave("./plots_and_tables/cpm_plots_top_DEGs_genrich.png", dpi = 300,
       width = 8, height = 4)

#### heatmaps
macs3_cpm_sig = cpm_macs3[c(macs3_sig$peak[macs3_sig$logFC > 0][1:75],
                            macs3_sig$peak[macs3_sig$logFC < 0][1:75]),]
row.names(macs3_cpm_sig) =paste(
  row.names(macs3_cpm_sig),
  macs3_sig$SYMBOL[match(row.names(macs3_cpm_sig), macs3_sig$peak)]
)

plt = macs3_cpm_sig %>%
  pheatmap(annotation_col = dge_macs3$samples %>%
             select(celltype,  drug),
           scale = "row",
           fontsize_row = 4,
           border_color = NA,
           treeheight_row = 0,
           
           treeheight_col = 20)
ggsave("./plots_and_tables/heatmap_DEGs_macs3.png", dpi = 300,
       plot = plt, width = 4, height = 8)


genrich_cpm_sig = cpm_genrich[c(genrich_sig$peak[genrich_sig$logFC > 0][1:75],
                                genrich_sig$peak[genrich_sig$logFC < 0][1:75]),]
row.names(genrich_cpm_sig) =paste(
  row.names(genrich_cpm_sig),
  genrich_sig$SYMBOL[match(row.names(genrich_cpm_sig), genrich_sig$peak)]
)

plt = genrich_cpm_sig %>%
  pheatmap(annotation_col = dge_genrich$samples %>%
             select(celltype,  drug),
           scale = "row",
           fontsize_row = 4,
           border_color = NA,
           treeheight_row = 0,
           
           treeheight_col = 20)
ggsave("./plots_and_tables/heatmap_DEGs_genrich.png", dpi = 300,
       plot = plt, width = 4, height = 8)

##### quick comparison
## tricky, because we want the logFC even for the genes that aren't sig,
## if they were sig in the other caller
macs3_stats = macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)
genrich_stats = genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)
macs3_stats %>% 
  group_by(SYMBOL) %>%
  summarize(macs3_logFC = mean(logFC),
            macs3_FDR = mean(FDR)) %>%
  select(macs3_logFC, macs3_FDR, SYMBOL) %>%
  full_join(genrich_stats %>% 
              group_by(SYMBOL) %>%
              summarize(genrich_logFC = mean(logFC),
                        genrich_FDR = mean(FDR)) %>%
              select(genrich_logFC, genrich_FDR, SYMBOL) ) %>%
  mutate(macs3_sig = abs(macs3_logFC) >1 & macs3_FDR < 0.05) %>%
  mutate(genrich_sig = abs(genrich_logFC) >1 & genrich_FDR < 0.05) %>%
  filter(macs3_sig | genrich_sig) %>%
  mutate(significant = case_when(
    macs3_sig & genrich_sig ~ "both",
    macs3_sig ~ "macs3 only",
    genrich_sig ~ "genrich only"
  )) %>%
  mutate(SYMBOL = ifelse(significant == "both", SYMBOL, "")) %>%
  ggplot(aes(x = macs3_logFC, y = genrich_logFC, label = SYMBOL, color = significant))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, color = "black", max.overlaps = 50, 
                           force = 10,segment.size = 0.2, segment.color = "gray")+
  labs(x = "macs3 log2FC", y = "genrich log2FC")+
  theme_bw(10)
ggsave("./plots_and_tables/scatterplot_sig_genes.png", dpi = 300,
       width = 6, height = 4)

## tables 
writexl::write_xlsx(list("macs3" = macs3_stats,
                         "genrich" = genrich_stats),
                    "./plots_and_tables/SR_vs_ST_significance_tables_and_annotations.xlsx")


####### 
# Now we do UK vs mock in the SR
#######

## subset the objects and re-calculate norm factors
dge_macs3 = dge_macs3[, metadata_sub$sample[metadata_sub$celltype == "renewing"]]
dge_genrich = dge_genrich[, metadata_sub$sample[metadata_sub$celltype == "renewing"]]

dge_macs3 = calcNormFactors(dge_macs3)
dge_genrich = calcNormFactors(dge_genrich)

design <- model.matrix(~ drug, data = dge_macs3$samples)

macs3_fit = glmQLFit(dge_macs3, design = design)
macs3_res = glmQLFTest(macs3_fit, coef = 2)
genrich_fit = glmQLFit(dge_genrich, design = design)
genrich_res = glmQLFTest(genrich_fit, coef = 2)


macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0) %>%
  pull(logFC_dir) %>%
  table() # 148 False, 29 TRue

genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0) %>%
  pull(logFC_dir) %>%
  table() # 99 false, 15 true



data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, "")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)+
  labs(x = "logFC (UK / mock)")
ggsave("./plots_and_tables/volcano_drug_effect_in_SR_macs3.png", dpi = 300,
       width = 6, height = 4)

library(tidyr)

data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(annotation = ifelse(FDR < 0.05 & abs(logFC) > 1,
                             annotation, NA)) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL,
             color = annotation))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)+
  labs(x = "logFC (UK / mock)")
ggsave("./plots_and_tables/volcano_drug_effect_in_SR_annotations_macs3.png", dpi = 300,
       width = 7, height = 4)



data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, "")) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)+
  labs(x = "logFC (UK / mock)")
ggsave("./plots_and_tables/volcano_drug_effect_in_SR_genrich.png", dpi = 300,
       width = 6, height = 4)


data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(annotation = ifelse(FDR < 0.05 & abs(logFC) > 1,
                             annotation, NA)) %>% 
  mutate(SYMBOL = ifelse(FDR < 0.05 & abs(logFC) > 1,
                         SYMBOL, NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(FDR), label = SYMBOL,
             color = annotation))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, max.overlaps = 20,
                           force = 10)+
  theme_bw(10)+
  labs(x = "logFC (UK / mock)")
ggsave("./plots_and_tables/volcano_drug_effect_in_SR_annotations_genrich.png", dpi = 300,
       width = 7, height = 4)

### bar charts of peak location annotations
data.frame(macs3_res) %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(abs(logFC) > 1 & FDR < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "higher in UK", "higher in mock")) %>%
  count(annotation, direction) %>%
  mutate(n = ifelse(direction == "higher in mock", -n, n)) %>%
  ggplot(aes(x = annotation, y = n, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("higher in mock" = "blue", "higher in UK" = "red")) +
  theme_bw(10) +
  labs(x = "Annotation", y = "# significant peaks", fill = "binding") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
ggsave("./plots_and_tables/barchart_drug_effect_in_SR_annotations_macs3.png", dpi = 300,
       width = 5, height = 2.5)



data.frame(genrich_res) %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  separate(annotation, into = c("annotation"),sep = "\\(") %>% 
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(abs(logFC) > 1 & FDR < 0.05) %>%
  mutate(direction = ifelse(logFC > 0, "higher in UK", "higher in mock")) %>%
  count(annotation, direction) %>%
  mutate(n = ifelse(direction == "higher in mock", -n, n)) %>%
  ggplot(aes(x = annotation, y = n, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("higher in mock" = "blue", "higher in UK" = "red")) +
  theme_bw(10) +
  labs(x = "Annotation", y = "# significant peaks", fill = "binding") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3)
ggsave("./plots_and_tables/barchart_drug_effect_in_SR_annotations_genrich.png", dpi = 300,
       width = 5, height = 2.5)


###### some point plots for top genes for gut check
cpm_macs3 = cpm(dge_macs3, log = TRUE)
cpm_genrich = cpm(dge_genrich, log = TRUE)
macs3_sig = macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0)  %>%  
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)

cpm_macs3[c(macs3_sig$peak[macs3_sig$logFC > 0][1:20],
            macs3_sig$peak[macs3_sig$logFC < 0][1:20]),] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "CPM") %>%
  left_join(metadata_sub) %>%
  left_join(macs3_sig) %>%
  separate(annotation, into = c("annotation"), sep = "\\(", extra = "drop") %>%
  arrange(logFC) %>%
  mutate(annotation = factor(annotation, levels = unique(annotation)),
         SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)),
         peak = factor(peak, levels = unique(peak))) %>%
  mutate(drug = ifelse(drug == "none", "mock", drug)) %>%
  ggplot(aes(x = drug, color = drug, y = CPM))+
  geom_point(size = 1.5,position = position_dodge(width = 0.2))+
  facet_wrap(peak ~ SYMBOL ~ annotation, scales = "free_y", nrow =4)+
  theme_bw(5.5)+
  scale_color_manual(values = c("mock" = "blue", "UK" = "red")) +
  labs(y = "log10(CPM)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(3, "mm"))
ggsave("./plots_and_tables/cpm_plots_drug_effect_in_SR_top_DEGs_macs3.png", dpi = 300,
       width = 8, height = 4)

genrich_sig = genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  mutate(logFC_dir = logFC > 0)  %>%  
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)

cpm_genrich[c(genrich_sig$peak[genrich_sig$logFC > 0][1:15],
              genrich_sig$peak[genrich_sig$logFC < 0][1:20]),] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "CPM") %>%
  left_join(metadata_sub) %>%
  left_join(genrich_sig) %>%
  separate(annotation, into = c("annotation"), sep = "\\(", extra = "drop") %>%
  arrange(logFC) %>%
  mutate(annotation = factor(annotation, levels = unique(annotation)),
         SYMBOL = factor(SYMBOL, levels = unique(SYMBOL)),
         peak = factor(peak, levels = unique(peak))) %>%
  mutate(drug = ifelse(drug == "none", "mock", drug)) %>%
  ggplot(aes(x = celltype, color = drug, y = CPM))+
  geom_point(size = 1.5,position = position_dodge(width = 0.2))+
  facet_wrap(peak ~ SYMBOL ~ annotation, scales = "free_y", ncol =10)+
  theme_bw(5.5)+
  scale_color_manual(values = c("mock" = "blue", "UK" = "red")) +
  labs(y = "log10(CPM)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(3, "mm"))
ggsave("./plots_and_tables/cpm_plots_drug_effect_in_SR_top_DEGs_genrich.png", dpi = 300,
       width = 8, height = 4)

#### heatmaps
macs3_cpm_sig = cpm_macs3[c(macs3_sig$peak[macs3_sig$logFC > 0][1:29],
                            macs3_sig$peak[macs3_sig$logFC < 0][1:75]),]
row.names(macs3_cpm_sig) =paste(
  row.names(macs3_cpm_sig),
  macs3_sig$SYMBOL[match(row.names(macs3_cpm_sig), macs3_sig$peak)]
)

plt = macs3_cpm_sig %>%
  pheatmap(annotation_col = dge_macs3$samples %>%
             mutate(drug = ifelse(drug == "none", "mock", "UK")) %>%
             select(celltype,  drug),
           scale = "row",
           fontsize_row = 4,
           border_color = NA,
           treeheight_row = 0,
           annotation_colors = list(drug = c("mock" = "blue", "UK" = "red")),
           treeheight_col = 20)
ggsave("./plots_and_tables/heatmap_drug_effect_in_SR_DEGs_macs3.png", dpi = 300,
       plot = plt, width = 4, height = 8)


genrich_cpm_sig = cpm_genrich[c(genrich_sig$peak[genrich_sig$logFC > 0][1:15],
                                genrich_sig$peak[genrich_sig$logFC < 0][1:75]),]
row.names(genrich_cpm_sig) =paste(
  row.names(genrich_cpm_sig),
  genrich_sig$SYMBOL[match(row.names(genrich_cpm_sig), genrich_sig$peak)]
)

plt = genrich_cpm_sig %>%
  pheatmap(annotation_col = dge_genrich$samples %>%
             mutate(drug = ifelse(drug == "none", "mock", "UK")) %>%
             select(celltype,  drug),
           scale = "row",
           fontsize_row = 4,
           border_color = NA,
           treeheight_row = 0,
           annotation_colors = list(drug = c("mock" = "blue", "UK" = "red")),
           treeheight_col = 20)
ggsave("./plots_and_tables/heatmap_drug_effect_in_SR_DEGs_genrich.png", dpi = 300,
       plot = plt, width = 4, height = 8)

##### quick comparison
## tricky, because we want the logFC even for the genes that aren't sig,
## if they were sig in the other caller
macs3_stats = macs3_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  rename(peak = genes) %>%
  left_join(macs3[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_macs3_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)
genrich_stats = genrich_res %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(PValue, method = "BH")) %>% 
  rename(peak = genes) %>%
  left_join(genrich[,1:5] %>%
              mutate(Chr = paste0("chr", Chr)) %>%
              rename(peak = Geneid,
                     seqnames = Chr,
                     start = Start, 
                     end = End)) %>%
  left_join(annotated_genrich_bed %>%
              data.frame()  
  ) %>%
  arrange(FDR)
macs3_stats %>% 
  group_by(SYMBOL) %>%
  summarize(macs3_logFC = mean(logFC),
            macs3_FDR = mean(FDR)) %>%
  select(macs3_logFC, macs3_FDR, SYMBOL) %>%
  full_join(genrich_stats %>% 
              group_by(SYMBOL) %>%
              summarize(genrich_logFC = mean(logFC),
                        genrich_FDR = mean(FDR)) %>%
              select(genrich_logFC, genrich_FDR, SYMBOL) ) %>%
  mutate(macs3_sig = abs(macs3_logFC) >1 & macs3_FDR < 0.05) %>%
  mutate(genrich_sig = abs(genrich_logFC) >1 & genrich_FDR < 0.05) %>%
  filter(macs3_sig | genrich_sig) %>%
  mutate(significant = case_when(
    macs3_sig & genrich_sig ~ "both",
    macs3_sig ~ "macs3 only",
    genrich_sig ~ "genrich only"
  )) %>%
  mutate(SYMBOL = ifelse(significant == "both", SYMBOL, "")) %>%
  ggplot(aes(x = macs3_logFC, y = genrich_logFC, label = SYMBOL, color = significant))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2, color = "black", max.overlaps = 50, 
                           force = 10,segment.size = 0.2, segment.color = "gray")+
  labs(x = "macs3 log2FC", y = "genrich log2FC")+
  theme_bw(10)
ggsave("./plots_and_tables/scatterplot__drug_effect_in_SR_sig_genes.png", dpi = 300,
       width = 6, height = 4)

## tables 
writexl::write_xlsx(list("macs3_drug_effect_in_SR" = macs3_stats,
                         "genrich_drug_effect_in_SR" = genrich_stats),
                    "./plots_and_tables/drug_effect_in_SR_significance_tables_and_annotations.xlsx")
