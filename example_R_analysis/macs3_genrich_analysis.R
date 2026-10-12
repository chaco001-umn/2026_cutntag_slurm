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


macs3 = read_table("/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/macs3/peak_counts", comment = "#")
genrich = read_table("/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/genrich/peak_counts", comment = "#")

# toss non-chromosomal
macs3 = macs3 %>%
  filter(!is.na(Chr))
names(macs3)[7:ncol(macs3)] = basename(names(macs3)[7:ncol(macs3)]) %>% 
  str_remove(".bam") %>%
  str_sub(1,6)

genrich = genrich %>%
  filter(!is.na(Chr))
names(genrich)[7:ncol(genrich)] = basename(names(genrich)[7:ncol(genrich)]) %>% 
  str_remove(".bam") %>%
  str_sub(1,6)

###### Compare the bed files

macs3_bed = macs3[,1:5] %>%
  tibble::column_to_rownames("Geneid")  %>%
  mutate(Chr = paste0("chr", Chr)) %>%
  GenomicRanges::makeGRangesFromDataFrame()
genrich_bed = genrich[,1:5] %>%
  tibble::column_to_rownames("Geneid")  %>%
  mutate(Chr = paste0("chr", Chr)) %>%
  GenomicRanges::makeGRangesFromDataFrame()



hits <- findOverlaps(macs3_bed, genrich_bed) # query, subject

ov <- pintersect(
  macs3_bed[queryHits(hits)],
  genrich_bed[subjectHits(hits)]
)

percent_of_macs3_peaks_overlapped_by_genrich <- width(ov) / width(macs3_bed[queryHits(hits)])
percent_of_genrich_peaks_overlapped_by_macs3 <- width(ov) / width(genrich_bed[subjectHits(hits)])
data.frame(metric = c("total overlapping peak pairs",
                      "macs3 peaks",
                      "genrich peaks",
                      "macs3 peaks also found in genrich",
                      "genrich peaks also found in macs3",
                      "median macs3 width",
                      "median genrich width",
                      "mean percent of macs3 peak covered by intersecting genrich peak",
                      "mean percent of genrich peak covered by intersecting macs3 peak"),
           value = c(nrow(data.frame(hits)),
                     nrow(macs3),
                     nrow(genrich),
                     length(unique(queryHits(hits))),
                     length(unique(subjectHits(hits))),
                     median(width(macs3_bed)),
                     median(width(genrich_bed)),
                     round(100 * mean(percent_of_macs3_peaks_overlapped_by_genrich)),
                     round(100 * mean(percent_of_genrich_peaks_overlapped_by_macs3)))
)  %>% 
  write_tsv("./plots_and_tables/caller_comparison_stats.tsv")

library(org.Hs.eg.db)
annotated_macs3_bed = ChIPseeker::annotatePeak(macs3_bed, 
                              TxDb =TxDb.Hsapiens.UCSC.hg38.knownGene,
                              annoDb = "org.Hs.eg.db")
annotated_genrich_bed = ChIPseeker::annotatePeak(genrich_bed, 
                               TxDb =TxDb.Hsapiens.UCSC.hg38.knownGene,
                               annoDb = "org.Hs.eg.db")
annotated_beds = list("macs3" = annotated_macs3_bed,
                      "genrich" = annotated_genrich_bed)
plotAnnoBar(annotated_beds)+
  theme_bw(10)+
  labs(title = "")+
  theme(legend.key.height = unit(2, "mm"))
ggsave("./plots_and_tables/plot_anno.png",dpi = 300, width = 6, height = 3)


plotDistToTSS(annotated_beds)+
  theme_bw(10)+
  labs(title = "")+
  theme(legend.key.height = unit(2, "mm"))
ggsave("./plots_and_tables/plot_disttotss.png",dpi = 300, width = 6, height = 3)

#### enrichment to get at types of genes
entrezids = data.frame(annotated_genrich_bed)$geneId %>% 
  as.character() %>% 
  unique()


ego <- enrichGO(gene = entrezids, 
                keyType = "ENTREZID", 
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
ego = simplify(ego)
dotplot(ego, showCategory = 50,font.size = 4, label_format = 100)
ggsave("./plots_and_tables/enrichment_of_genrich_peaks.png",dpi = 300, width = 6, height = 6)



###### Now do stats

metadata = read_tsv("/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/samples.txt", col_names = c("sample", "trash", "treatment",
                                                   "trash2", "date")) %>%
  as.data.frame() %>%
  mutate(treatment = str_remove(treatment, "CT29 ")) %>% 
  select(sample, treatment, date) %>%
  mutate(celltype = ifelse(str_detect(treatment, "SR"),
                           "renewing", "synciotroph")) %>%
  mutate(drug = ifelse(str_detect(treatment, "UK"),
                       "UK", "none")) 

dge_macs3 = DGEList(macs3 %>%
                tibble::column_to_rownames("Geneid") %>%
                {.[, metadata$sample]},genes = macs3$Geneid,
              samples = metadata %>%
                tibble::column_to_rownames("sample"))


dge_genrich = DGEList(genrich %>%
                      tibble::column_to_rownames("Geneid") %>%
                      {.[, metadata$sample]},genes = genrich$Geneid,
                    samples = metadata %>%
                      tibble::column_to_rownames("sample"))

### filter low "expression"
keep_macs3 <- filterByExpr(dge_macs3, group = dge_macs3$samples$celltype,min.count = 20)
keep_genrich <- filterByExpr(dge_genrich, group = dge_genrich$samples$celltype,min.count = 20)
table(keep_macs3)
table(keep_genrich)

### what does intersection look after filtering?
macs3_sub = macs3[keep_macs3,]
genrich_sub = genrich[keep_genrich,]


macs3_bed_sub = macs3_sub[,1:5] %>%
  tibble::column_to_rownames("Geneid")  %>%
  mutate(Chr = paste0("chr", Chr)) %>%
  GenomicRanges::makeGRangesFromDataFrame()
genrich_bed_sub = genrich_sub[,1:5] %>%
  tibble::column_to_rownames("Geneid")  %>%
  mutate(Chr = paste0("chr", Chr)) %>%
  GenomicRanges::makeGRangesFromDataFrame()

hits <- findOverlaps(macs3_bed_sub, genrich_bed_sub) # query, subject
findOverlaps(macs3_bed, genrich_bed) 


#### now start some unsupervised stuff
dge_macs3 <- dge_macs3[keep_macs3, , keep.lib.sizes = FALSE]
macs3 = macs3[keep_macs3, ]

dge_genrich <- dge_genrich[keep_genrich, , keep.lib.sizes = FALSE]
genrich = genrich[keep_genrich, ]



dge_macs3 = calcNormFactors(dge_macs3)
dge_genrich = calcNormFactors(dge_genrich)

plt = cpm(dge_macs3) %>%
  cor %>%
  {.^2} %>%
  pheatmap(annotation_col = dge_macs3$samples %>%
             select(celltype, drug, date),
           border_color = NA,
           treeheight_row = 0,
           treeheight_col = 20)
ggsave("./plots_and_tables/corrmap_macs3.png",
       plot = plt,dpi =300,
       width = 5, height = 4)

plt = cpm(dge_genrich) %>%
  cor %>%
  {.^2} %>%
  pheatmap(annotation_col = dge_genrich$samples %>%
             select(celltype, drug, date),
           border_color = NA,
           treeheight_row = 0,
           treeheight_col = 20)
ggsave("./plots_and_tables/corrmap_genrich.png",
       plot = plt,dpi =300,
       width = 5, height = 4)

cpmy_macs3 = cpm(dge_macs3,log = TRUE)
cpmy_genrich = cpm(dge_genrich,log = TRUE)
cpmy_variance_genrich = apply(cpmy_genrich, 1, sd)
cpmy_variance_macs3 = apply(cpmy_macs3, 1, sd)


plt = cpmy_macs3[cpmy_variance_macs3 >= sort(cpmy_variance_macs3, decreasing = TRUE)[500], ]%>%
  pheatmap(annotation_col = dge_macs3$samples%>%
             select(celltype, drug, date),
           show_rownames = FALSE, 
           scale = "row",
           border_color = NA,
           treeheight_row = 0,
           treeheight_col = 20)
ggsave("./plots_and_tables/vargenes_heatmap_macs3.png",
       plot = plt,dpi =300,
       width = 4, height = 4)


plt = cpmy_genrich[cpmy_variance_genrich >= sort(cpmy_variance_genrich, decreasing = TRUE)[500], ]%>%
  pheatmap(annotation_col = dge_genrich$samples %>%
             select(celltype, drug, date),
           show_rownames = FALSE, 
           scale = "row",
           border_color = NA,
           treeheight_row = 0,
           treeheight_col = 20)
ggsave("./plots_and_tables/vargenes_heatmap_genrich.png",
       plot = plt,dpi =300,
       width = 4, height = 4)


######pca
macs3_pca = cpmy_macs3[cpmy_variance_macs3 >= sort(cpmy_variance_macs3, decreasing = TRUE)[5000], ]%>%
  t() %>%
  prcomp(scale. = TRUE,
         center = TRUE) 
macs3_pca$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, label = sample,color = celltype, shape = drug))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2)+
  labs(x = paste("PC1: ",round(macs3_pca$sdev[1]), "%"),
       y = paste("PC2: ",round(macs3_pca$sdev[2]), "%"),)+
  theme_bw(10)
ggsave("./plots_and_tables/pca_macs3.png",dpi =300,
       width = 6, height = 4)

genrich_pca = cpmy_genrich[cpmy_variance_genrich >= sort(cpmy_variance_genrich, decreasing = TRUE)[5000], ]%>%
  t() %>%
  prcomp(scale. = TRUE,
         center = TRUE) 
genrich_pca$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, label = sample,color = celltype, shape = drug))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2)+
  labs(x = paste("PC1: ",round(genrich_pca$sdev[1]), "%"),
       y = paste("PC2: ",round(genrich_pca$sdev[2]), "%"),)+
  theme_bw(10)
ggsave("./plots_and_tables/pca_genrich.png",dpi =300,
       width = 6, height = 4)

### pcs 2,3
macs3_pca$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC2, y = PC3, label = sample,color = celltype, shape = drug))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2)+
  labs(x = paste("PC2: ",round(macs3_pca$sdev[2]), "%"),
       y = paste("PC3: ",round(macs3_pca$sdev[3]), "%"),)+
  theme_bw(10)
ggsave("./plots_and_tables/pca23_macs3.png",dpi =300,
       width = 6, height = 4)


genrich_pca$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC2, y = PC3, label = sample,color = celltype, shape = drug))+
  geom_point()+
  ggrepel::geom_text_repel(size = 2)+
  labs(x = paste("PC2: ",round(genrich_pca$sdev[2]), "%"),
       y = paste("PC3: ",round(genrich_pca$sdev[3]), "%"),)+
  theme_bw(10)
ggsave("./plots_and_tables/pca23_genrich.png",dpi =300,
       width = 6, height = 4)


#### go through all that again but skip the non 04/30 samples
metadata_sub = metadata %>%
  filter(date == "04/30/24")
macs3 = read_table("/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/macs3/peak_counts", comment = "#")
genrich = read_table("/projects/standard/swernimo/shared/RI/2025_cutntag/Project_002_output/genrich/peak_counts", comment = "#")

macs3 = macs3 %>%
  filter(!is.na(Chr))
names(macs3)[7:ncol(macs3)] = basename(names(macs3)[7:ncol(macs3)]) %>% 
  str_remove(".bam") %>%
  str_sub(1,6)

genrich = genrich %>%
  filter(!is.na(Chr))
names(genrich)[7:ncol(genrich)] = basename(names(genrich)[7:ncol(genrich)]) %>% 
  str_remove(".bam") %>%
  str_sub(1,6)
dge_macs3 = DGEList(macs3 %>%
                      tibble::column_to_rownames("Geneid") %>%
                      {.[, metadata_sub$sample]},genes = macs3$Geneid,
                    samples = metadata_sub %>%
                      tibble::column_to_rownames("sample"))


dge_genrich = DGEList(genrich %>%
                        tibble::column_to_rownames("Geneid") %>%
                        {.[, metadata_sub$sample]},genes = genrich$Geneid,
                      samples = metadata_sub %>%
                        tibble::column_to_rownames("sample"))

### filter low "expression"
keep_macs3 <- filterByExpr(dge_macs3, group = dge_macs3$samples$celltype,min.count = 20)
keep_genrich <- filterByExpr(dge_genrich, group = dge_genrich$samples$celltype,min.count = 20)
table(keep_macs3)
table(keep_genrich)

dge_macs3 = dge_macs3[keep_macs3, , keep.lib.sizes = FALSE]
dge_genrich = dge_genrich[keep_genrich, , keep.lib.sizes = FALSE]

dge_macs3 = calcNormFactors(dge_macs3)
dge_genrich = calcNormFactors(dge_genrich)


save(annotated_macs3_bed, annotated_genrich_bed,
     genrich, genrich_bed,
     macs3, macs3_bed,
     dge_genrich, dge_macs3,metadata,metadata_sub,
     file = "./analysis/intermediate_dataobjects.Rdat")
rm(list = ls())
##########
# pause point
#
#########
