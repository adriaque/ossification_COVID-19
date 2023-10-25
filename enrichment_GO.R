setwd('/Users/alejandroadriaquelozano/Documents/Systems Biology/Network Biology/Project/Modules/')

library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)

#exporting pre-processed data for lung and colon
module_1 <- read_csv('Module_0_1 default node.csv')
module_2 <- read_csv('Module_0_2 default node.csv')
module_4 <- read_csv('Module_0_4 default node.csv')


test_1 <- enrichGO(module_1$`shared name`, keyType = 'SYMBOL' ,ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)
test_2 <- enrichGO(module_2$`shared name`, keyType = 'SYMBOL' ,ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)
test_4 <- enrichGO(module_4$`shared name`, keyType = 'SYMBOL' ,ont = 'BP', OrgDb = 'org.Hs.eg.db', pvalueCutoff = 1, qvalueCutoff = 1)

summary_1 <- as.data.frame(summary(test_1))
summary_1$GeneRatio<-paste0(" ",summary_1$GeneRatio)
write_excel_csv(summary_1,'results_GO_1.csv')
dotplot(test_1, showCategory =10)

summary_2 <- as.data.frame(summary(test_2))
summary_2$GeneRatio<-paste0(" ",summary_2$GeneRatio)
write_excel_csv(summary_2,'results_GO_2.csv')
dotplot(test_2, showCategory =10)

summary_4 <- as.data.frame(summary(test_4))
summary_4$GeneRatio<-paste0(" ",summary_4$GeneRatio)
write_excel_csv(summary_4,'results_GO_4.csv')
dotplot(test_4, showCategory =10)
