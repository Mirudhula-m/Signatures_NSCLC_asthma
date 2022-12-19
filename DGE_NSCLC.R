# Script to get data from GEO

library(GEOquery)
library(dplyr) # querying columns in subject data
library(limma)
library(ggplot2)
library(ggrepel)
library(readr)

###################
# NSCLC Dataset
###################

GEOaccession <- 'GSE29013'
gse <- getGEO(GEOaccession)
gse <- gse[[1]]

boxplot(exprs(gse),outline=FALSE)


# Get the subject information of interest
sampleInfo <- pData(gse)
# making a separate group column
groupCol <- data.frame(group=c('Cancer'))
for(i in 1:(nrow(sampleInfo)-1)) { groupCol <- rbind(groupCol, list('Cancer')) }
sampleInfo <- bind_cols(sampleInfo, groupCol)
sampleInfo <- select(sampleInfo, group)

# heatmap of the gene expression data
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)  


###################
# Normal Dataset
###################

# Get normal subjects dataset, shared between Asthma dataset
# Normal/Asthma data
GEOaccession <- 'GSE64913'
gse_norm <- getGEO(GEOaccession)
gse_norm <- gse_norm[[1]]

# Get the subject information of interest
sampleInfo_norm <- pData(gse_norm)
sampleInfo_norm <- select(sampleInfo_norm, `diagnosis:ch1`)
sampleInfo_norm <- rename(sampleInfo_norm, group=`diagnosis:ch1`)

# heatmap of the gene expression data
corMatrix <- cor(exprs(gse_norm),use="c")
pheatmap(corMatrix, annotation_col=sampleInfo_norm)  

# Filter out patients who are healthy in the asthma dataset
filtered_patient_ID = row.names(subset(sampleInfo_norm,group=='Healthy'))
filtered_sampleInfo_norm <- subset(sampleInfo_norm,group=='Healthy')
filtered_expr <- exprs(gse_norm)[,c(filtered_patient_ID)]

# heatmap of the gene expression data
corMatrix <- cor(filtered_expr,use="c")
pheatmap(corMatrix)

# Combined data
combined_expr <- bind_cols(exprs(gse), filtered_expr) 
combined_sampleInfo <- bind_rows(sampleInfo, filtered_sampleInfo_norm) 
rownames(combined_expr) <- rownames(fData(gse))
  write_csv(combined_expr, "combined_expression_lung.csv")
corMatrix <- cor(combined_expr,use="c")
#rownames(combined_sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix, annotation_col=combined_sampleInfo)     
boxplot(combined_expr,outline=FALSE)

# DGE of NSCLC against Normal patients
# Create one-hot encoding
design <- model.matrix(~0+combined_sampleInfo$group)

# Renaming column names
colnames(design) <- c("Tumour","Normal")

## calculate median expression level
cutoff <- median(as.numeric(unlist(combined_expr)))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- combined_expr > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
# combined_expr_sub <- combined_expr[keep,]
g1 <- gse
filtered_gse <- g1[keep,]
filtered_expr_only_lc <- exprs(g1)
selected_IDs <- row.names(filtered_expr_only_lc)
combined_expr_sub <- combined_expr %>% filter(row.names(combined_expr) %in% c(selected_IDs))

fit <- lmFit(combined_expr_sub, design)
head(fit$coefficients)

# defining the contrast
contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)

# Apply empirical Bayes step to perform DGE
fit2 <- eBayes(fit2)

table(decideTests(fit2))

# annotations
anno <- fData(gse)
anno <- select(anno,ID,`Gene Symbol`,ENTREZ_GENE_ID)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"Gene Symbol")

# cutoffs
p_cutoff <- 10e-60
fc_cutoff <- 5
topN <- 5

# identify dge
# add a column of NAs
full_results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
full_results$diffexpressed[full_results$logFC > fc_cutoff & full_results$adj.P.Val < p_cutoff] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
full_results$diffexpressed[full_results$logFC < -fc_cutoff & full_results$adj.P.Val < p_cutoff] <- "DOWN"


# plot dge
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, FC = abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=diffexpressed))  + geom_point() + scale_color_manual(values=c("#FF9999","#000033", "#00CC66"))


# Filtering the data for the top 20
# Get the upregulated genes set as a csv
up = filter(full_results, adj.P.Val < p_cutoff, logFC > fc_cutoff)
fin_up = data.frame(up$ID)
write_csv(fin_up, file="./LC_geneset_up.csv")

# Get the downregulated genes set as a csv
down = filter(full_results, adj.P.Val < p_cutoff, logFC < -fc_cutoff)
fin_down = data.frame(down$ID)
write_csv(fin_down, file="./LC_geneset_down.csv")


# Get the gene set as a csv
filter(full_results, adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff) %>%
  write_csv(file="filtered_de_results_new.csv")


t =topTable(fit2, number=32977)
# t is the list of all the genes with p values and log fc
write_csv(t, file="./LC_all_genes_p_logfc.csv")

