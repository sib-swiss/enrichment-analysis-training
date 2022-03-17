## **Bonus code** :champagne_glass:

The following code was added thanks to questions from course participants of past sessions. They might be useful for you too.


## Gene label conversions


It is often useful to convert different types of gene labels even further. Here is an example for gene label conversion 
and gene information extraction using [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html). A lot of information for each gene can be obtained, such as chromosome location,
description, biotype, or the symbols of mouse homologs of human genes, etc.

```r
# install and load the package:
BiocManager::install("biomaRt")
library(biomaRt)
# list the available options
listEnsembl()

ensembl <- useEnsembl(biomart = "genes")

# List the available species and reference genomes:
datasets <- listDatasets(ensembl)
head(datasets)
# search for human dataset and genome version:
datasets[grep("sapiens", datasets$dataset),]
# 80 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

# To obtain information on human genes, first the data specific to human has to be accessed from the online Ensembl repository
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# you can extract gene biotype, description, wiki gene description,
# chromosome location of genes, even strand information etc
# using the listAttributes function allows you to view the type of information per gene you can extract.
attributes <- listAttributes(ensembl)
attributes[1:15,]

# Now we extract symbol (external_gene_name), description, gene_biotype (eg whether it is a protein coding or other type of gene):
biomart_gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", "wikigene_description"),
                           filter="ensembl_gene_id",
                           values=NK_vs_Th$ensembl_gene_id,
                           mart=ensembl)
# check the structure of the resulting data frame:
head(biomart_gene_info)

# search for the oncogene key word in the description:
biomart_gene_info[grep("oncogene", biomart_gene_info$description),]
biomart_gene_info[grep("tumor protein", biomart_gene_info$description),]

# search for genes that have TP53 in their gene symbol:
biomart_gene_info[grep("TP53", biomart_gene_info$external_gene_name),]
```

For conversion of human gene symbols to mouse homologs for example (or vice versa if you provide mouse genes as the "values" arguments
                                                                    and "hsapiens_homolog_associated_gene_name" in the list of the attributes argument), you can also use biomaRt.

```r
# Convert human genes to mouse homologs
ensembl_human_to_mouse <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "mmusculus_homolog_associated_gene_name"),
                                filter="ensembl_gene_id",
                                values=NK_vs_Th$ensembl_gene_id,
                                mart=ensembl)
head(ensembl_human_to_mouse)
# ensembl_gene_id external_gene_name mmusculus_homolog_associated_gene_name
# 1 ENSG00000000003             TSPAN6                                 Tspan6
# 2 ENSG00000000419               DPM1                                   Dpm1
# 3 ENSG00000000419               DPM1                                Gm20716
# 4 ENSG00000000457              SCYL3                                  Scyl3
# 5 ENSG00000000460           C1orf112                               BC055324
# 6 ENSG00000000938                FGR                                    Fgr
```

## msigdbr package

The [msigdbr](https://cran.r-project.org/web/packages/msigdbr/index.html) package hosted on CRAN allows to access gene set collections hosted on MSigDB directly within R. Check out its vignette to view how to download collections for other species such as mouse, etc. The function ```msigdbr_species()``` allows you to list available species. 

For example, to download the Hallmark collection with human gene symbols within R:
```r
gmt <- msigdbr::msigdbr(species = "human", category = "H")

# Create the 2 column-format (TERM2GENE argument) required by clusterProfiler:
h_gmt <- gmt[,c("gs_name", "gene_symbol")]

```


## Code for a heatmap of p-values

For heatmaps, ggplot2 can also be used. Here is an example for a heatmap of the p-values of 6 different gene sets (gs), and the
p-value for each gene set was calculated in two comparisons, so we compared the enrichment in genes differentially expressed
between NK cells and Th cells, and between NK cells and CD8 T cells (note that this is a dummy example just to show you the 
                                                                     example of the code, it is not based on real RNA seq data).
Before using ggplot2, you need to create a dataframe that contains a column with the gene set identities, a column with the 
name of the cell type comparisons, and a column with the p-value of each gene set in each comparison. So a specific format is 
required for ggplot2.

```r
# Create a data frame that contains the p-value for every gene set for every 
# cell type comparison. You need to include the values also for the non-significant
# gene sets. If you use function of the clusterProfiler package, you will usually
# obtain results in the @result slot only for significant gene sets. But if you need
# p-values also for the non-significant gene sets, you can change the argument
# pvalueCutoff = 1 in the gseGO function (and other functions of the clusterProfiler
# package)

# For the example, create a dummy data frame with the list of gene sets (gs), the
# 2 cell type comparisons, and the p-value for each gene set in each comparison:
ora_to_plot<-as.data.frame(cbind(ID=c("gs1", "gs2","gs3","gs4", "gs5", "gs6",
                                      "gs1", "gs2","gs3","gs4", "gs5", "gs6"),
                                 comparison=c(rep("NK_vs_Th",6), rep("NK_vs_CD8",6)),
                                 p_val=as.numeric(c(0.8,0.9,0.6, 0.054, 0.00001, 0.0003,
                                                    0.0002, 0.004, 0.001,0.01, 0.85,0.9))))

# transform the p-val to -log10(p-val):
ora_to_plot$log10_p_val<--log10(as.numeric(ora_to_plot$p_val))
# set the none-significant p-values to 0 so that they appear grey in the heatmap:
ora_to_plot$log10_p_val<-ifelse(ora_to_plot$log10_p_val<(-log10(0.05)), 0, ora_to_plot$log10_p_val)
range(ora_to_plot$log10_p_val) # [1]  0 5

head(ora_to_plot, n=12)
#     ID comparison p_val log10_p_val
# 1  gs1   NK_vs_Th   0.8    0.000000
# 2  gs2   NK_vs_Th   0.9    0.000000
# 3  gs3   NK_vs_Th   0.6    0.000000
# 4  gs4   NK_vs_Th 0.054    0.000000
# 5  gs5   NK_vs_Th 1e-05    5.000000
# 6  gs6   NK_vs_Th 3e-04    3.522879
# 7  gs1  NK_vs_CD8 2e-04    3.698970
# 8  gs2  NK_vs_CD8 0.004    2.397940
# 9  gs3  NK_vs_CD8 0.001    3.000000
# 10 gs4  NK_vs_CD8  0.01    2.000000
# 11 gs5  NK_vs_CD8  0.85    0.000000
# 12 gs6  NK_vs_CD8   0.9    0.000000

# create a palette of colors based on the Plasma palette (grDevices package)
plot(1:20, 1:20, col=hcl.colors(20,"Plasma"), pch=15, cex=3)
hcl.colors(9,"Plasma")# "RdYlGn")
breaks<-seq(from=0, to=5, by=0.5 )

color<-c("grey78", rev(hcl.colors(9,"Plasma")))

# create plot and export as png:
p<-ggplot(ora_to_plot, aes(comparison, ID, fill= log10_p_val)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(breaks= breaks, 
                       colors = color) +
  theme_bw() 
ggsave(plot = p, filename = "heatmap_p_value_ORA.png", 
       device="png",
       width = 4,height = 6)
```

You will obtain the following heatmap:
  
  <figure>
  <img src="../assets/images/heatmap_p_value_ORA.png" width="300"/>
  </figure>