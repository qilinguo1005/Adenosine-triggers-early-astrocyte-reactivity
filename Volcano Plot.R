
data<-read.csv('D:/0_GLACxA1AR/Figure for paper/RNA-Seq/ip/6h.Hisat2.DESeq2.cko_vs_ctl.csv')

data<-na.omit(data)
data<-data[!duplicated(data$gene_name),]


logFC_t <- with(data,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange))) 


logFC_t <- round(logFC_t, 3)
logFC_t =0.3


data$Change = as.factor(ifelse(data$padj< 0.05 & abs(data$log2FoldChange) > logFC_t,
                               
                               ifelse(data$log2FoldChange > logFC_t ,'UP','DOWN'),'STABLE'))
table(data$Change) 

this_tile <- paste0('Volcano plot cKO vs ctl 6h',
                    '\nCutoff for Log2FC is ',round(logFC_t,3),
                    '\nThe number of up genes is ',nrow(data[data$Change =='UP',]) ,
                    '\nThe number of down genes is ',nrow(data[data$Change =='DOWN',]))

install.packages("rlang")
library(ggplot2)
ggplot(data, aes(x=log2FoldChange, y=-log10(padj),color=Change)) + 
  geom_point(alpha=0.4, size=2) +  # 设置点的透明度和大小
  theme_bw(base_size = 12) +  #设置一个主题背景
  xlab("Log2(Fold change)") + # x轴名字
  ylab("-Log10(P.adj)") + # y轴名字
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('steelblue','gray','brown')) + # 各自的颜色
  geom_hline(yintercept = -log10(0.05), lty = 4) + #定义p值和线形
  geom_vline(xintercept = c(-logFC_t, logFC_t), lty = 4)+ #定义差异倍数和线形
  labs(title = this_tile) #加上题目

data$label <- ifelse(data$padj< 0.05& abs(data$log2FoldChange) >= 1,data$gene_name,"")

data$label <- ifelse(data$gene_name %in% b$gene_name, data$gene_name, "")


library(ggrepel)
ggplot(data, aes(x=log2FoldChange, y=-log10(padj),color=Change)) + 
  geom_point(alpha=0.4, size=2) + 
  theme_bw(base_size = 12) + 
  xlab("Log2(Fold change)") +
  ylab("-Log10(P.adj)") +
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','gray50','red')) +
  geom_hline(yintercept = -log10(0.05), lty = 4) +
  geom_vline(xintercept = c(-logFC_t, logFC_t), lty = 4)+
  labs(title = this_tile)+
  xlim(-5.5,4.5)+
  ylim(0,10)+
  geom_text_repel(data = data, aes(label = label),size=4.5,max.overlaps = 1000)+
  theme_classic(base_size = 16)
  geom_label_repel(data = data, aes(label = label),
                   size = 3,box.padding = unit(0.25, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 1000)
