
#################################################
####                                         ####
####      cluster profiler                   ####           
####                                         ####
#################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DOSE)

## GO enrichment 
### dot GO:

setwd("F:/TCGA_yasaman/")

y <- read.delim("F:/TCGA_yasaman/for cluster.txt", sep="\t",header = T)

y2 <- y$gene_name
ent_uni <- bitr(y2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ent_uni <- ent_uni$ENTREZID
ego_BP <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "BP" ,
                   #universe = ent_uni 
)

ego_MF <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "MF" ,
                   #universe = ent_uni 
)

ego_CC <- enrichGO(gene = ent_uni,
                   OrgDb = org.Hs.eg.db ,
                   ont = "CC" ,
                   #universe = ent_uni 
)



#N <- as.numeric(sub("\\d+/", "", ego[1, "BgRatio"]))


p1 <- dotplot(ego_BP, showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Biological Process") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))



p2 <- dotplot(ego_MF, showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Molecular Function") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))


p3 <- dotplot(ego_CC,  showCategory = 10
              #,x = ~Count/N
) +
  #ggplot2::xlab("Rich Factor")+
  ggtitle("Cellular Component") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+xlim(0.01, 0.09)+
  guides("Corrected P-Value" = guide_legend(order = 1),
         size = guide_legend(order = 2,title = "Gene number"))+
  theme(plot.title = element_text(hjust=0.5,size = 15, face = "bold.italic"))


png("dotplot_yasaman.png",width = 2500, height = 4000,res = 300,units = "px")


cowplot::plot_grid(p1, p2, p3, nrow =3,
                   labels=c("I","II","III"))

dev.off()


### cnet GO:

geneList <- y$logFC
names(geneList) <- as.character(y$symbol)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)


edox1 <- setReadable(ego_BP, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)

edox2 <- setReadable(ego_MF, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)

edox3 <- setReadable(ego_CC, 'org.Hs.eg.db', 'ENTREZID')
view(edox1)


k1 <- cnetplot(edox1,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Biological Process") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


k2 <- cnetplot(edox2,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Molecular Function") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


k3 <- cnetplot(edox3,colorEdge = TRUE,categorySize=2,
               foldChange = geneList) +
  ggtitle("Cellular Component") +
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold.italic"))


png("cnetplot_GO.png",width = 10000, height = 3000,res = 300,units = "px")

cowplot::plot_grid(k1, k2, k3, ncol =3,
                   labels=c("I","II","III"))

dev.off()

####################################################



## KEGG enrichment :

ekg <- enrichKEGG(ent_uni,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

#deactive pvalue and qvalue effect on the numbers of pathways in results (get all results)

view(ekg)

#view gene symbols:

ekgm<- setReadable(ekg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#cALCULATE RICH FACTOR:
ekgm <- mutate(ekgm, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#make data frame to read by ggplot2:
ekgm <- as.data.frame(ekgm)
#write enrichKEGG
write.table(ekgm, file="enrichKEGG.txt", quote=F, sep="\t")


####### dotplot BY Rich factor:

# select first 20:
ekgm <-arrange(ekgm, p.adjust) %>% slice(1:20)

#tiff("dotplot_KEGG.tiff",width = 3000, height = 2000,res = 300,units = "px")

d2 <- ggplot(ekgm,
             aes(richFactor,Description)) + 
  #geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE),name="Corrected P-Value") +
  scale_colour_gradientn(
    colours = rainbow(10, s = 1, v = 1, start = 0, end = max(1, 10 - 1)/10,),
    name="Corrected P-Value")+
  scale_size_continuous(range=c(2, 6)) +
  #theme_minimal() + 
  theme_bw(base_size = 14)+ 
  xlab("Rich factor") +
  ylab(NULL) +  guides(size=guide_legend("Gene number"))+
  ggtitle("Statistics of KEGG Enrichment")

#dev.off()

####################################################

### cnet KEGG:
geneList <- x$logFC
names(geneList) <- as.character(x$geneSymbol)
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)


jj <-arrange(ekg@result, p.adjust)
ekg@result <- jj

ekgx <- setReadable(ekg, 'org.Hs.eg.db', 'ENTREZID')
view(ekgx)


#tiff("cnetplot_KEGG.tiff",width = 4300, height = 3500,res = 300,units = "px")

p2 <- cnetplot(ekgx,colorEdge = TRUE,categorySize=2,
               foldChange = geneList,circular = TRUE) + 
  scale_color_gradient2(name="logFc",low="#0C8B21", mid="white",
                        high="#E50303", space ="Lab" )+
  guides("size" = guide_legend(order = 1),
         "logFC" = guide_legend(order = 2),
         "category" = guide_legend(order = 3))+
  theme(plot.title = element_text(hjust=0.5,size = 30, face = "bold"))


#dev.off()

#COMBINE CNET:
tiff("cnetplot_combine.tiff",width = 8500, height = 3500,res = 300,units = "px")

plot_grid(
  p2 ,NULL,p1,
  rel_widths = c(1, 0.05, 1),
  nrow = 1,
  labels = c('A','', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_colour = "black",
  label_size = 30
)

dev.off()

#COMBINE DOTPLOT:
tiff("dotplot_combine.tiff",width = 5500, height = 2300,res = 300,units = "px")

plot_grid(
  d2 ,d1,
  rel_widths = c(1.75, 2),
  nrow = 1,
  labels = c('A', 'B'), 
  label_fontfamily = "serif",
  label_fontface = "plain",
  label_colour = "black",
  label_size = 30
)

dev.off()
##############################

### Enrichment Map using KEGG :
ekg <- enrichKEGG(ent_gene,
                  pvalueCutoff = 1)



tiff("EnrichmentMap_KEGG.tiff",width = 4000, height = 4000,res = 300,units = "px")

emapplot(ekg,showCategory = 30,pie="count", pie_scale=1.5, layout="kk")

dev.off()

