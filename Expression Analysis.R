
#################################################
####                                         ####
####     Gene Expression Data Analysis       ####
####                                         ####
#################################################

setwd("Your Working Directory")

#### Determining tumoral and non-tumoral samples

exp.N <- colData(exp)$definition %in% c("Solid Tissue Normal")
exp.T <- colData(exp)$definition %in% c("Primary solid Tumor")


#### Pre-processing and normalizing expression data

exp.prep <- TCGAanalyze_Preprocessing(object = exp, cor.cut = 0.6)

exp.norm <- TCGAanalyze_Normalization(tabDF = exp.prep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

exp.filt <- TCGAanalyze_Filtering(tabDF = exp.norm,
                                  method = "quantile",
                                  qnt.cut = 0.25)

#### Extraction of differentially-expressed genes

exp.DEGs <- TCGAanalyze_DEA(mat1 = exp.filt[,exp.T],
                            mat2 = exp.filt[,exp.N],
                            Cond1type = "Primary solid Tumor",
                            Cond2type = "Solid Tissue Normal",
                            pipeline = "edgeR",
                            method = "exactTest",
                            fdr.cut = 1,
                            logFC.cut = 0)

write.table(exp.DEGs, "Differentially-expressed Genes.txt", quote = F)

#### Drawing the volcano plot

TCGAVisualize_volcano(x = exp.DEGs$logFC,
                      y = exp.DEGs$FDR,
                      x.cut = 1,
                      y.cut = 0.05,
                      width = 15,
                      height = 10,
                      legend = "State",
                      color = c("black","red","green"),
                      xlab = "Gene expression fold change (Log2)",
                      title = "Volcano plot (Primary solid Tumor vs Solid Tissue Normal)",
                      filename = "Volcano plot of STAD-DEGs.pdf",
                      show.names = F) 

