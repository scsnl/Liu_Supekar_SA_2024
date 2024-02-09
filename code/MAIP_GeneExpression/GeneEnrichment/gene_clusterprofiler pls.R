# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# install.packages("org.Hs.eg.db", repos="http://bioconductor.org/packages/3.2/data/annotation")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

d = read.csv(file="~/Documents/projects/2021_Longt_math_gene/final/outputs/GeneExpression/PLS1_geneWeights_descend.csv",header = FALSE)
geneList=d[,3]
names(geneList)=as.character(d[,1])
geneList = sort(geneList,decreasing = TRUE)
geneSet = read.gmt("math_Lee_new1.gmt")
res<-GSEA(geneList,exponent = 1,
          eps = 1e-10,
          #nPerm = 10000,
          minGSSize = 3,
          maxGSSize = 750,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          TERM2GENE=geneSet,
          TERM2NAME = NA,
          verbose = TRUE,
          seed = FALSE,
          by = "fgsea"
)
summary(res)
GeneSetOrder<-res$ID
NES_value<-res$NES
Padjust<-res$p.adjust
Results<-data.frame(GeneSetOrder,NES_value,Padjust)

save(res, file="res.RData")
save(Results,file="Results.RData")
write.csv(Results, file = "Results.csv",row.names = FALSE)


######################################################################

# visualizing the results of GSEA

#load("res.RData")

gseaplot2(res, 
          geneSetID="MathAbility",
          title="Math ability", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('MathAbility.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

gseaplot2(res, 
          geneSetID="HighestMath",
          title="Highest math class taken", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('HighestMath.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

gseaplot2(res, 
          geneSetID="Wordreading2020",
          title="Wordreading 2020", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Wordreading2020.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

gseaplot2(res, 
          geneSetID="Wordreading2022",
          title="Wordreading 2022", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Wordreading2022.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

gseaplot2(res, 
          geneSetID="Dyslexia2022",
          title="Dyslexia 2022", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('dyslexia2022.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

gseaplot2(res, 
          geneSetID="WM",
          title="Working memory", 
          color = "#2F5C85",
          base_size = 25,
          rel_heights = c(0.8, 0.2, 0.5),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")
ggsave('Working memory2014.tiff', width =15,height = 12,dpi = 300, limitsize = FALSE,units = "cm")

