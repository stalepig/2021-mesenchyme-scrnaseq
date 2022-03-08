setwd("/Volumes/LEBRON/Chicago/03_10_21_Safina_mesenchyme_scRNA_Seq/")

library(monocle3)
library(ggplot2)
source("/Volumes/SHAQ/2019/01_15_19_DSS_single_cell_post/monocle3/camMonocleHelper.R")

# cds.samp1 <- load_cellranger_data(pipestance_path = "BP10/",genome = "mm10")
# cds.samp2 <- load_cellranger_data(pipestance_path = "BP11/",genome = "mm10")
# cds.samp3 <- load_cellranger_data(pipestance_path = "BP12/",genome = "mm10")
# cds.samp4 <- load_cellranger_data(pipestance_path = "BP13/",genome = "mm10")
# cds.samp5 <- load_cellranger_data(pipestance_path = "BP14/",genome = "mm10")
# cds.samp6 <- load_cellranger_data(pipestance_path = "BP15/",genome = "mm10")

# cds.full <- combine_cds(cds_list = list(cds.samp1,cds.samp2,cds.samp3,cds.samp4,cds.samp5,cds.samp6))

# cds.full <- recluster.cds(cds = cds.full)

# mtr <- top_markers(cds = cds.full,group_cells_by = "cluster",genes_to_test_per_group = 25,verbose = T)
# goi <- filter.marker.results(marker_test_res = mtr,num.to.keep = 4,criterion = "specificity")
# g <- plot_genes_by_group(cds = cds.full,markers = goi,max.size = 5)
# print(g)

# g10 <- plot_cells(cds = cds.full,show_trajectory_graph = F,group_label_size = 6)
# print(g10)

# g20 <- plot_cells(cds = cds.full,show_trajectory_graph = F,label_cell_groups = F,genes = c("Tnfrsf1a","Tnfrsf1b","Tnf","Bmp4","Vim","Wnt5a","Pdgfra","Acta2","Col6a1","Fn1","Gata6","Birc5","Ccnb1","Krt8","Gli1","Rspo3"))
# print(g20)

# colData(cds.full)$genotype <- ifelse(test = colData(cds.full)$sample %in% c(1,2,3),yes = "wt",no = "ko")
# colData(cds.full)$cell.group <- clusters(cds.full)
# colData(cds.full)$umap.1 <- reducedDims(x = cds.full)[["UMAP"]][,1]
# colData(cds.full)$umap.2 <- reducedDims(x = cds.full)[["UMAP"]][,2]

# g30 <- plot_cells(cds.full,show_trajectory_graph = F,label_cell_groups = F,color_cells_by = "genotype")
# print(g30)

# barcodes <- colData(cds.full) %>% data.frame
# taba <- table(barcodes$sample,barcodes$cell.group)

## generates assignments
# vec.clus <- c("1-ISEMF",
#                  "2-ISEMF",
#                  "3-ISEMF",
#                  "4-Generic",
#                  "5-ISEMF",
#                  "6-ISEMF",
#                  "7-ISEMF",
#                  "8-Generic",
#                  "9-Generic",
#                  "10-ISEMF,Mmp9",
#                  "11-Generic",
#                  "12-ISEMF",
#                  "13-Generic",
#                  "14-Generic",
#                  "15-Fragments")
# assignments <- data.frame(Key=vec.clus)
# assignments$Cluster <- sapply(X = as.character(assignments$Key),FUN = function(str) strsplit(x = str,split = "-")[[1]][1])
# assignments$Assignment <- sapply(X = as.character(assignments$Key),FUN = function(str) strsplit(x = str,split = "-")[[1]][2])
# cds.full <- assign.clusters(cds = cds.full,assignments = assignments,column = 3,new.col.name = "annotation")
colData(cds.full)$condition <- paste(colData(cds.full)$annotation,colData(cds.full)$genotype,sep = "-")

## finds marker genes for pathway analysis - enrichr
# mtra <- top_markers(cds = cds.full,group_cells_by = "condition",genes_to_test_per_group = 200,verbose = T)

## removes fragments genes because the cluster is too small to make reliable calls
# mtra10 <- subset(mtra,!grepl(pattern = "Fragments",x = mtra$cell_group,ignore.case = T))

## makes the enrichr table
# mtra20 <- subset(mtra10,mtra10$specificity>0.30)
# mtra25 <- mtra20[order(mtra20$specificity,decreasing = T),]
# mtra30 <- subset(mtra25,!duplicated(mtra20$gene_short_name))
# write.table(x = mtra30,file = "table_markers.csv",sep = ",",row.names = F)

## makes a gsea table
# barcodes <- colData(cds.full) %>% data.frame
# cds.full.filtered <- cds.full[,subset(barcodes,barcodes$annotation!="Fragments") %>% row.names]
# agg <- aggregate.expression.by.factor(cds = cds.full.filtered,grouping = "condition",do.print = T)
# d <- data.frame(agg)
# d$id <- row.names(d)
# d10 <- melt(data = d,id.vars = "id")
# colnames(d10) <- c("id","cell_group","abundance")
# feats <- fData(x = cds.full.filtered) %>% data.frame
# d20 <- merge(x = d10,y = feats,by="id",all = F)
# d30 <- aggregate(formula = abundance~gene_short_name+cell_group,data = d20,FUN = sum)
# d40 <- dcast(data = d30,formula = gene_short_name~cell_group,value.var = "abundance")
# hgene <- read.csv(file = "homologene.csv",header = F)
# colnames(hgene) <- c("ID","Organism","V3","Gene_symbol","V5","V6")
# hgene$Gene_symbol <- as.character(hgene$Gene_symbol)
# hgene10 <- dcast(data = hgene,formula = ID~Organism,value.var = "Gene_symbol",fun.aggregate = function(x) {x[1]})
# hgene20 <- hgene10[,c("ID","9606","10090")]
# colnames(hgene20) <- c("ID","Human","Mouse")
# d50 <- merge(x = hgene20,y = d40,by.x = "Mouse",by.y = "gene_short_name",all.x = F,all.y = T)
# d60 <- na.omit(d50)
# d70 <- d60[,-c(1)]
# d80 <- d70[,c(2,1,3:7)]
# colnames(d80)[1:2] <- c("Name","Description")
# write.table(x = d80,file = "gsea_full_filtered.txt",sep = "\t",row.names = F,quote = F)