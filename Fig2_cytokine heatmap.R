library(pheatmap)
data<-read.csv("Fig2_cpa_as.csv",head = T, row.names=1)
pheatmap(log2(data), scale="none",
         color = colorRampPalette(c("#08519C", "white", "#A50F15"))(16),
         border="white", 
         cluster_cols = F, 
         cluster_rows = F,
         show_rownames = T,
         fontsize_row = 10,
         fontsize_col = 20,
         angle_col = 90,
         breaks = seq(-2,2,by = 0.25)
)
