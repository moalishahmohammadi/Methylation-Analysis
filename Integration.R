#################################################
####                                         ####
####      integration and visualization      ####           
####                                         ####
#################################################

setwd("Your Directory")
#adj.p,symbol,logFc,deltaBeta,symbol

deg<- read.delim("F:/TCGA_yasaman/diff with genes name.txt")


dmg<- read.delim("F:/TCGA_yasaman/data_common.txt")


deg<-deg[,c(2,5,6)]

dmg<-dmg[,c(2,3,4)]

dmg1=aggregate(x =dmg,by=list(dmg$Gene_Symbol), max )
deg1=aggregate(x =deg,by=list(deg$Gene_Symbol), max )


# names in dmg are different . they have more space so delete the space by below cod

dmg$Gene_Symbol <- gsub(' ','',dmg$Gene_Symbol)

data<-merge(deg,dmg,by.x= "Gene_Symbol" , by.y= "Gene_Symbol")


data$Group <- ifelse(data$logFC <(-1) & data$`delta.beta` <(-0.2) & 
                       data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hypo-Down" ,
                     ifelse(data$logFC >1 & data$`delta.beta` <(-0.2) &
                              data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hypo-Up",
                            ifelse(data$logFC <(-1) & data$`delta.beta` >0.2 &
                                     data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hyper-Down",
                                   ifelse(data$logFC >1 & data$`delta.beta` >0.2 &
                                            data$FDR<0.05 & data$`adj.pvalue`<0.05, "Hyper-Up",
                                          "Not-sig"))))


#write.csv(data,"final.csv")
#x<-data[data$Group !="Not-sig",]

write.table(data,"data.txt",
            quote = F,sep = "\t")
# plotting

library(ggplot2)

cols <- c("Hypo-Down" = "yellow", "Hypo-Up" = "blue", "Not-sig" = "grey", "Hyper-Down" = "red", "Hyper-Up" = "springgreen4")

ggplot(data, aes(x = `delta.beta`, y =logFC , color = Group)) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = cols) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_hline(yintercept = -1, linetype="dashed") + 
  geom_vline(xintercept = 0.2, linetype="dashed") + 
  geom_vline(xintercept = -0.2,  linetype="dashed") +
  xlab("Mean methylation diffrences") + 
  ylab("Log2 expression change") 

