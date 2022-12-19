
library(GEOquery)
library(dplyr) # querying columns in subject data
library(limma)
library(ggplot2)
library(ggrepel)
library(readr)

###################
# Asthma Dataset
###################

GEOaccession <- 'GSE64913'
gse <- getGEO(GEOaccession)
gse <- gse[[1]]

pData(gse) ## print the sample information
P <- pData(gse)
P_healthy = subset(P, `diagnosis:ch1`=="Healthy")
P_asthma = subset(P, `diagnosis:ch1`=="Severe Asthmatic")
fData(gse) ## print the gene annotation

boxplot(exprs(gse),outline=FALSE)

#Inspecting clinical variables
library(dplyr)
sampleInfo <- pData(gse)

## renaming to more convenient column names
sampleInfo <- rename(sampleInfo, group = `diagnosis:ch1`, patient = `patient id:ch1`)

#Sample clustering and PCA
library(pheatmap)
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix) 

BiocManager::install("pheatmap")

## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix, annotation_col=sampleInfo) 

#DGE
library(limma)
design <- model.matrix(~0+sampleInfo$group)

#Renaming the columns
colnames(design) <- c("Healthy","Asthma")

#Filtering lowly expressed genes
summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)


contrasts <- makeContrasts(`Healthy` - `Asthma`, levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

topTable(fit2)
decideTests(fit2)
table(decideTests(fit2))

aw <- arrayWeights(exprs(gse),design)
fit <- lmFit(exprs(gse), design,weights = aw)
contrasts <- makeContrasts(`Healthy` - `Asthma`, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
decideTests(fit2)
table(decideTests(fit2))

anno <- fData(gse)
anno <- select(anno,'ID',`Gene Symbol`,`ENTREZ_GENE_ID`)
anno <- rename(anno, Gene = `Gene Symbol`)
fit2$genes <- anno
t =topTable(fit2, number=32977)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 0.75

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
  ggplot(aes(x = logFC, y = B, col=diffexpressed)) + geom_point() + scale_color_manual(values=c("#FF9999","#000033", "#00CC66"))

# Get the upregulated genes set as a csv
up = filter(full_results, adj.P.Val < p_cutoff, logFC > fc_cutoff)
fin_up = data.frame(up$ID)
write_csv(fin_up, file="./Asthma_geneset_up.csv")

# Get the downregulated genes set as a csv
down = filter(full_results, adj.P.Val < p_cutoff, logFC < -fc_cutoff)
fin_down = data.frame(down$ID)
write_csv(fin_down, file="./Asthma_geneset_down.csv")

# t is the list of all the genes with p values and log fc
write_csv(t, file="./Asthma_all_genes_p_logfc.csv")



