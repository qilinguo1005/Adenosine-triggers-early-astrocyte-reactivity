library(pheatmap)

select1<-read.csv("C:/Users/Qilin/Desktop/select1.csv",head = T, row.names=1)


data<-log2(select1+1)

p<- pheatmap(data, scale="row",
         color = colorRampPalette(c("blue", "white", "red"))(256),
         border="white", 
         cluster_cols = F, 
         cluster_rows = T,
         show_rownames = T,
         cutree_rows =6
         )

