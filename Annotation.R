
#################################################
####                                         ####
####      Annotation of expression data      ####        
####                                         ####
#################################################

genes<- read.delim("gene-report.txt")

genes<-substr(genes$gene_id, 1 , 15 )

genes=aggregate(x =genes,by=list(genes$x), max )

row.names(genes)<- genes$x


exp.DEGs$gene_name<- genes[as.character(rownames(exp.DEGs)),"gene_name"]

exp.DEGs$gene_type<- genes[as.character(rownames(exp.DEGs)),"gene_type"]

write.table(exp.DEGs,"diff with genes name.txt",
            quote = F,sep = "\t")




#################################################
####                                         ####
####      Annotation of methylation data     ####      
####                                         ####
#################################################

annot.met <- as.data.frame(read.delim("annotation.txt"))

annot.met<-annot.met[,c("Composite.Element.REF","Gene_Symbol")]

dmc<- as.data.frame(met.state[,c("mean.Primary.solid.Tumor.minus.mean.Solid.Tissue.Normal",
                                 "p.value.adj.Primary.solid.Tumor.Solid.Tissue.Normal")])

colnames(dmc)<-c("delta beta","adj pvalue")



data_common <- merge(dmc,annot.met,by.x="row.names",by.y="Composite.Element.REF")


data_common <- tidyr::separate_rows(data_common, Row.names, Gene_Symbol,sep = ";") 

data_common$comb <- paste(data_common$Row.names,data_common$Gene_Symbol)

data_common <- data_common[!duplicated(data_common$comb), ]

data_common <- data_common[, -5]

write.table(data_common,"data_common.txt",
            quote = F,sep = "\t")