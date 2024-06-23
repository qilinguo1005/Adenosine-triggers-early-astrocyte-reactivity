###GSEA

install.packages("devtools")   # unnecessary if you have it already
devtools::install_github("Bioconductor/BiocManager", ref="ghost-binary-repo")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

###Install and load required packages
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("GSEABase")
install.packages("GSEABase", dependencies = TRUE)

library(GSEABase)
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

#SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

###Prepare Input
# reading in data from deseq2
df = read.csv("cKO_vs_ctl.csv",header = TRUE)
df = df[1:13]
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$gene_name

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

###Gene Set Enrichment
#ont one of “BP”, “MF”, “CC” or “ALL”
#nPerm the higher the number of permutations you set, the more accurate your result will, but the longer the analysis will take.
#minGSSize minimum number of genes in set (gene sets with lower than this many genes in your dataset will be ignored).
#maxGSSize maximum number of genes in set (gene sets with greater than this many genes in your dataset will be ignored).
#pvalueCutoff pvalue Cutoff.
#ENSEMBL SYMBOL
#pAdjustMethod one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 20, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr")

selected_pathways <- gse$Description
selected_pathways

dotplot(gse, showCategory=30,label_format=100,font.size = 8,
        title = "GSEA-GO term analysis: ip_6h cKO vs ctl", split=".sign") + facet_grid(.~.sign,scales = "free")+scale_size(range=c(1.5, 5))
p2 <- p + scale_color_continuous(low='red', high='blue')
p2
dotplot(gse, showCategory = selected_pathways[y])

### Encrichment Map
x2 <-enrichplot::pairwise_termsim(gse)
x2 <- pairwise_termsim(gse) 
x2$Description
yy= c(212,272)

emapplot(x2, showCategory = pathways)
emapplot(gse, showCategory = selected_pathways[y])

### Category Net #plotcategorySize can be either 'pvalue' or 'geneNum'
y= c(6,7,9,10,18,27)
y= c(135,119,283,229,54,423)
cnetplot(gse, categorySize="geneNum", foldChange=gene_list, showCategory = selected_pathways[y])


###GSEA Plot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 6)


###KEGG Gene Set Enrichment Analysis
##Prepare Input
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"), OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$gene_name %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = T)

write.csv(kegg_gene_list,"C:/Users/Qilin/Desktop/kegg_gene_list_24h.csv")

###Create gseKEGG object
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 20,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

selected_pathways2 <- kk2$Description
selected_pathways2
z= c(3,4,6,7,8,9,10,15,18,19,20,23,26,27,28,30,32,33)


###Dotplot
dotplot(kk2, showCategory = selected_pathways2[z] , label_format=100,  title = "Enriched Pathways" ,orderBy = "x", font.size = 10,decreasing = TRUE, split=".sign",) + facet_grid(.~.sign)+scale_size(range=c(1.5, 5))
