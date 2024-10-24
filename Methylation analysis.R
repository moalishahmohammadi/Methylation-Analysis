
#################################################
####                                         ####
####        Methylation data analysis        ####
####                                         ####
#################################################
setwd("your Directory")
#### Pre-processing 
#### Removing NAs, sex chromosomes and SNPs

met <- subset(met, subset = (rowSums(is.na(assay(met))) == 0))
met <- subset(met, subset = !as.character(seqnames(met))
              %in% c("*", "chrX", "chrY"))

#### Extraction of differentially-methylated CpGs.

met.state <- TCGAanalyze_DMC(met,
                             groupCol = "definition",
                             group1 = "Primary solid Tumor",
                             group2 = "Solid Tissue Normal",
                             p.cut = 0.05,
                             diffmean.cut = 0.2,
                             legend = "State",
                             plot.filename = "Volcano plot of STAD-DMCs.pdf")

met.DMCs <- subset(met.state,
                   subset = met.state$status != "Not Significant")

met.DMC.hypo <- subset(met.state,
                       met.state$status == "Hypomethylated in Primary solid Tumor")

met.DMC.hyper <- subset(met.state,
                        met.state$status == "Hypermethylated in Primary solid Tumor")
