
################################################################### #
#                                                                   #
#    版本信息：clustering and DR(通用)    发布时间：2020-03-30      #
#    Pipeline versiion: 2.10                                        #
#                                                                   #
################################################################### # 



#===========================0. 载入需要的工具包======================================================

rm(list = ls())#清除内存变量

library(flowCore)
library(Rcpp)
library(cytofkit)
library(igraph)
library(ggplot2)
library(ggthemes)
library(Rtsne)
library(dplyr)
library(cytofexplorer)
library(Rmisc)
library(stringi)
library(RColorBrewer)

#===========================1.Pipeline设置==========================================================

#聚类算法：

    Do_PhenoGraph=TRUE
    Do_FlowSOM=TRUE

#降维算法：

    Do_tSNE=TRUE
    DO_UMAP=TRUE

#===========================2. 目录设置==============================================================

#  (A) Input：<--
      raw_fcs_dir="01_rawfcs"  #设置读取数据的文件夹名称

#  (B) 自动设置文件夹
      wdir<-dirname(rstudioapi::getActiveDocumentContext()$path) #自动获取projectdir
      projectdir<-sub("/[^/]*$","",wdir)
      raw_fcs_dir<-paste0(projectdir,"/",raw_fcs_dir)
      metadata_dir<-paste0(projectdir,"/02_metadata")
      setwd(wdir)


#  (C)checkpoint：(!注意：目录中不要有中文字符)
      wdir
      projectdir
      source(paste0(projectdir,"/03_backups/Version_control.R"))
      
#===========================3. 数据读取（读取->Downsample->合并）=======================================

#  (A) Input：<--
      mergeMethod <- "ceil" #合并方法："ceil", "all", "fixed", "min"
      fixedNum <- 5000       # 设置从每个文件抽取的细胞数目


#  (B) 读取文件，DownSample，Merge
      file_name <- list.files(raw_fcs_dir,pattern='.fcs$', full=TRUE)
      combined_data_raw <- cytof_exprsMerge(fcsFiles = file_name,
                                      transformMethod = "none",
                                      mergeMethod =mergeMethod,
                                      fixedNum = fixedNum)
#  (C) 简化列名
      paraname_raw<-colnames(combined_data_raw)
      paraname<-sub("<NA>","",paraname_raw)
      paraname<-sub(".*<","",paraname)
      paraname<-sub(">.*","",paraname)
      paraname<-sub("-","_",paraname)
      paraname<-sub("^[0-9][0-9][0-9][A-Z][a-z]_|^89Y_","",paraname)  #去掉数字标签
      paraname<-sub("_EQ[0-9]$","",paraname)  #去掉EQ标签

      colnames(combined_data_raw)<-paraname

#  (D) 增加File_ID
      File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
      combined_data_raw<-data.frame(combined_data_raw,File_ID)

#  (E) Checkpoint
      head(combined_data_raw)

#===========================4. 追加MetaData数据========================================================

#  (A) marker选择，仅仅需要input一次
#     1) 生成marker选择文件
    #  markers2csv(paraname_raw,colnames(combined_data_raw))

#     2) Input：<--
#      到metadata目录中找到all_markers.csv,
#      --检查markers一列名称是否OK，如果需要，可以直接修改；
#      --在各列上标记出需要选择的通道（标记1），保存；
#      --(*各列设置可以在后面进行，也可以根据需要添加新的列，注意列名要与相应的功能模块匹配。)

#     3) 重新读取保存的csv文件：
      all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),stringsAsFactors =FALSE)
      colnames(combined_data_raw)<-all_markers[,"markers"]

#  
#  (B) 文件分组，仅仅需要input一次
#     1) 生成分组文件
      groups2csv(unique(File_ID))

#     2) Input：<--
#      手工设置分组信息（仅首次需要）
#      到metadata目录中找到samplename.csv,打开后在B列Short_name，输入每个样本的简称;
#      样本间Shortname不能有重复，且尽量简洁以便出现在输出的图片中；
#      C列以后根据实验设计，输入实际的分组情况例如：Timepiont, Tissue_Type,Patient_ID ect.

#     3) 重新读取保存的csv文件：
      groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors=FALSE)
      groups$File_ID<-unique(File_ID)

#  (C)checkpoint: -->  查看groups内容,确认分组信息被正常读取
      all_markers
      groups  


#===========================5. Transform数据转换========================================================

#  (A)Input: <--
#     （仅需设置一次，如前面设置过，直接运行此模块代码）
#     通道选择，all_markers.csv,在transform列上标记出要进行转换的通道（标记1），保存


#  (B) 数据转换
      all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
      transform_id=which(all_markers$transform==1)
      combined_data_transformed<-combined_data_raw
      combined_data_transformed[, transform_id] <- apply(combined_data_transformed[, transform_id,drop = FALSE],
                                                         2,cytofAsinh)
#  (C) checkpoint
      head(combined_data_transformed)

      
#===========================6. 运行 PhenoGraph聚类=====================================================

if(Do_PhenoGraph){      
      
#  (A) Input: <--  PhenoGraph参数设置：
#    1) 设定k值

  k=30  #计算Knn网络的“邻居”个数

#    2) 通道选择：到mall_markers.csv,在要进行PhenoGraph聚类的标签行标记1，保存。

#  (B) 数据准备
  all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
  PhenoGraph_id=which(all_markers$PhenoGraph==1)
  PhenoGraph_input_data=combined_data_transformed[,PhenoGraph_id]

#  (C) PhenoGraph聚类

#PG_elbow(PhenoGraph_input_data) #<--去掉行首#可以进行elbow测试,耗时较长，建议在total event较少的情况下使用

  PhenoGraph_result <-as.numeric(membership(Rphenograph(data = PhenoGraph_input_data,k=k)))
#  (D) Checkpoint:
  hist(PhenoGraph_result,unique(PhenoGraph_result))

}

    save(PhenoGraph_result,file="PhenoGraph_result")
  load("PhenoGraph_result")    

#===========================7.运行FlowSOM聚类===============================================================
if(Do_FlowSOM){
  
library(FlowSOM)
#  (A).Input： 
#    1) 设定xdim和ydim的数值
#      FloSOM第一步SOM聚类，需要预先指定cluster数目，最终会生成xdim*ydim个cluster，
xdim=20     #<-- SOM聚类会生成xdim*ydim个cluster
ydim=20     #<-- 建议xdim和ydim 取数值相等（或相近）的自然数

#    2)在all_marker中新增一列，命名为FlowSOM，选中的marker用1标记。

#  (B).运行SOM聚类：

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
colnames_raw<-colnames(combined_data_raw)
all_markers$markers<-colnames_raw[1:nrow(all_markers)]
FlowSOM_id=which(all_markers$FlowSOM==1)
FlowSOM_marker<-all_markers[FlowSOM_id,"markers"]
FlowSOM_input_data=combined_data_transformed[,FlowSOM_id]
head(FlowSOM_input_data)

map <- SOM(data= as.matrix(FlowSOM_input_data),
           xdim=xdim,
           ydim=ydim,
           silent = F)

SOM_result<-map$mapping[,1]

#  (C).运行metacluster
#     FlSOM第二步是对获得的SOM cluster再进行聚类，产生metacluster，数目可人为设定或者通过elbow_test得到

FlowSOM_combined<-data.frame(FlowSOM=SOM_result,
                             FlowSOM_input_data)

metacluster_result<-metaClustering(FlowSOM_combined,
                                   clustername = "FlowSOM",
                                   metaClustering_method = "metaClustering_PhenoGraph", #<-- 聚类方法
                                   k_value=10, #<-- 聚类方法中的k值
                                   elbow_test=T, #<-- 决定是否进行elbow test，当kvalue=NULL时，会自行选择进行elbowtest
                                   seed=123)     #<-- 设定seed可以保持tsne图保持一致

#   (D).Checkpoint
head(metacluster_result)

}
save(metacluster_result,file="metacluster_result")      
load(file="metacluster_result")     
      
#===========================7a.（*可选步骤）:设置cluster别名（Setting alias of clusters）=============================================
#  (A).生成cluster_alias csv： 
       generate_cluster_alias(clustering_result= PhenoGraph_result,
                       cluster_name="PhenoGraph")

#  (B).进入output文件夹找到生成的csv文件，手动设置第二列的列名和各个cluster的别名；      
      
#===========================7b.（*可选步骤）:手动合并亚群(Manually merge clusters)=======================

    if(TRUE){    #<--如果要进行合并括号中改为TRUE
    
      merge_result<-cluster_merge(cluster_result=PhenoGraph_result, #<--待合并的聚类结果  
                    merge_list=list(c(4,8,13),#Treg
                                    c(10,18),#cd4tcm
                                    c(15,22),#cd4tn
                                    c(16,20),#DNT
                                    c(2,23),#nkt
                                    c(14,21),#cd8tcm
                                    c(7,12,24),#cd8tn
                                    c(3,17,19)#tcrgd
                                    ),               #<--指出要合并的亚群例如：  list(c(1,2,3),c(6:9))
                    cluster_rename=T)                               #<--是否对cluster重新命名
      
      #Checkpoint
      hist(merge_result,unique(merge_result))
    }

save(merge_result,file="merge_result")

#===========================8.（*可选步骤） Cluster-marker表达Preview(Density和Heatmap)=============================
##整合数据，从各个Meta cluster或文件中抽取定量细胞进行可视化(Downsample)

#  (A) Input: 整合数据：
      combined_rawdata_analysed <- data.frame(combined_data_raw,
                                        metacluster=metacluster_result) #<-- 加入待整合的聚类结果
      
#  (B) 生成Marker-Cluster Heatmap和Density图片

      cluster_marker_preview(combined_rawdata_analysed,
                             groups = groups,
                             all_markers=all_markers,
                             cluster_color=dif_seq_rainbow,  #选择颜色 其他选择：rainbow，dif_seq_rainbow,brewer_color_sets
                             cluster_name="metacluster",      #cluster 通道的名称
                             colored_by="expression"         #"expression"
                             #cluster_id=NULL                # 选择要显示的cluster
                             )
      
      
     head( combined_rawdata_analysed)
      
      
#===========================9. 将结果导出成FCS文件=====================================================
      
#  (A) Input： 整合待导出数据

      combined_data_output <- combined_data_raw[,0]
      if(Do_PhenoGraph)  combined_data_output <- data.frame(combined_data_output,PhenoGraph = PhenoGraph_result)  #<-- 输入待导出的通道
      if(Do_FlowSOM)     combined_data_output <- data.frame(combined_data_output,metacluster = metacluster_result)  #<-- 输入待导出的通道
      
      
#  (B) CheckPoint:
      
      head(combined_data_output)

#  (C) 根据File_ID将合并数据还原成单个样本数据，导出FCS文件
      cytof_addToFCS_modified(
        data=combined_data_output,
        rawFCSdir=raw_fcs_dir,
        analyzedFCSdir=paste(wdir,"/output",sep = ""),
        newHeader = "mc_" #<-- 可以自行设定输出文件的前缀
        )

#===========================10. DownSample==============================================================  

#  (A) Input： 整合rawdata和聚类结果
      combined_rawdata_analysed <- combined_data_raw
      if(Do_PhenoGraph)  combined_rawdata_analysed <- data.frame(combined_rawdata_analysed,PhenoGraph = PhenoGraph_result)  #<-- 输入待导出的通道
      if(Do_FlowSOM)     combined_rawdata_analysed <- data.frame(combined_rawdata_analysed,metacluster = metacluster_result)  #<-- 输入待导出的通道
      
      

#  (B) downsample         
      combined_rawdata_sampled<-equal_sample(x=combined_rawdata_analysed,
                                          sample_method="ceil",         #<-- 两种取值："ceil", "all"
                                          sample_type="by_cond",        #<-- "by_cluster","by_file","by_cond"
                                          sample_num=50000,              #<-- 每个cluster/file/cond抽取的细胞数
                                          cluster_name="PhenoGraph",    #<-- Cluster的名称
                                          groups=groups,                #仅在by_cond时需要，其他时候可以注释掉
                                          sample_cond = "Group_1"   #仅在by_cond时需要，其他时候可以注释掉
                                          )
      combined_transdata_sampled<-combined_rawdata_sampled
      combined_transdata_sampled[, transform_id] <- apply(combined_transdata_sampled[, transform_id,drop = FALSE],
                                                         2,cytofAsinh)
#  (C) Checkpoint
      head(combined_transdata_sampled)
      head(combined_rawdata_sampled)

      
      
#===========================11a. 运行tSNE(BH-sne)降维(本步骤需要时间较长)==================================================
if(Do_tSNE){
      
#  (A) Input: <-- 
#    1) tSNE参数设置：
      max_iter=2500  #迭代次数
      perplexity=30  #困惑度
      seed=1500      #随机数种子
      theta=0.5      #权衡速度与准确度，越小越精确，越大速度越快
      dims = 2       #降维输出维度（默认是2）

#    2) 通道选择：到mall_markers.csv,找到tSNE一列，在要进行tSNE降维的marker标记1，保存。

#  (B) tsne分析:

      if (exists('seed')) set.seed(seed)
      all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
      tSNE_para_id=which(all_markers$tSNE==1)
      tSNE_input_data=combined_transdata_sampled[,tSNE_para_id]
      tsne_result <- Rtsne(tSNE_input_data,
                           initial_dims = ncol(tSNE_input_data),
                           pca = FALSE,
                           dims = dims,
                           check_duplicates = FALSE,
                           perplexity=perplexity,
                           max_iter=max_iter,
                           theta=theta)$Y
      row.names(tsne_result)<-row.names(tSNE_input_data)
      colnames(tsne_result)<-c("tsne_1","tsne_2")

#   (C) checkpoint: 绘制出降维分析的草图
      
      plot(tsne_result)
}   
      save(tsne_result,file="tsne_result")
      load("tsne_result")
      
     
#===========================11b. 降维结果可视化========================================================
if(Do_tSNE){
        
#  (A) Input：
#    1) 将rawdata和降维结果整合：
      combined_rawdata_plot   <- data.frame(combined_rawdata_sampled,
                                      tsne_result)           #<-- 降维结果
#  (B) Checkpoint:
      head(combined_rawdata_plot)
      #head(combined_transdata_plot)
      #combined_data_plot$PhenoGraph = factor(combined_data_plot$PhenoGraph,level=c(1,2,3,4,5,6,7,8,9,10,11),
      save(combined_rawdata_plot,file="combined_rawdata_plot")                                      # labels = c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11"))    
      
#  (C) 生成tsne-Phenograph系列图片
      draw_tsne_figs(combined_rawdata_plot,  
                     groups=groups,
                     cluster_color=dif_seq_rainbow,
                     cluster_name="PhenoGraph",
                     alias_as_label = F,
                     major_cond="Group_1",
                     cluster_labels_size=15,
                     dot_size = 2,
                     show_cluster_labels = T,
                     reduction_dm1="tsne_1",
                     reduction_dm2="tsne_2",
                     show_contour=F,
                     output_format = "pdf")

#  (D)生成tsne热图
#    1) Input：<--
      #通道选择：打开mall_markers.csv,找到/创建tsne_heatmap一列，在在要生成热图的所在行标记1，保存退出

#    2)生成tsne-heatmap：
      all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
      all_markers$markers<-colnames(combined_data_raw)[1:nrow(all_markers)]
      heatmap_tsne_id<-which(all_markers$tsne_heatmap==1)
      heatmap_tsne_markers=colnames(combined_data_transformed[,heatmap_tsne_id])
      
      draw_tsne_heatmaps(combined_rawdata_plot,  #注意：该函数使用的是rawdata
                         heatmap_tsne_markers=heatmap_tsne_markers,
                         trans_method="simpleAsinh",
                         groups=groups,
                         #major_cond="Tissue_Type",
                         dot_size=4,
                         #dot_colours = colorRampPalette(c("#4575B4" ,"#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59" ,"#D73027")),
                         show_legend=T,
                         show_scale=F,
                         show_contour = F,
                         output_format="tiff",
                         single_file=F)

}   
#===========================12a.运行UMAP降维分析=====================================================
if(DO_UMAP){
      
#rUMAP：
#install.packages("uwot")
#Input: UMAP参数设置：
library(uwot)
set.seed(123)
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
UMAP_id=which(all_markers$UMAP==1)
UMAP_input_data=combined_transdata_sampled[,UMAP_id]
UMAP_input_data=as.matrix(UMAP_input_data)

rumap_result <- umap(UMAP_input_data, n_neighbors = 15, learning_rate =1, min_dist=0.4,init = "random",n_epochs =200)

colnames(rumap_result)<-c("UMAP1","UMAP2")
head(rumap_result)
plot(rumap_result)

}

#===========================12b. 降维结果可视化========================================================
if(DO_UMAP){
#  (A) Input：
#    1) 将rawdata和降维结果整合：
combined_rawdata_plot   <- data.frame(combined_rawdata_sampled,
                                      rumap_result)           #<-- 降维结果
#  (B) Checkpoint:
head(combined_rawdata_plot)

#  (C) 生成tsne-Phenograph系列图片
draw_tsne_figs(combined_rawdata_plot,  
               groups=groups,
               cluster_color=dif_seq_rainbow,
               cluster_name="PhenoGraph",
               alias_as_label = F,
               major_cond="Group_1",
               cluster_barplot_data=metacluster_result,
               cluster_labels_size=15,
               dot_size = 3,
               show_cluster_labels = T,
               reduction_dm1="UMAP1",
               reduction_dm2="UMAP2",
               output_dir="cluster_UMAP_plots",
               #group_seq=c("PBMC","Biopsy"),
               show_contour=F,
               output_format = "pdf")


#  (D)生成tsne热图
#    1) Input：<--
#通道选择：打开mall_markers.csv,找到/创建tsne_heatmap一列，在在要生成热图的所在行标记1，保存退出

#    2)生成tsne-heatmap：
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
all_markers$markers<-colnames(combined_data_raw)[1:nrow(all_markers)]
heatmap_tsne_id<-which(all_markers$tsne_heatmap==1)
heatmap_tsne_markers=colnames(combined_data_transformed[,heatmap_tsne_id])

draw_tsne_heatmaps(combined_rawdata_plot,  #注意：该函数使用的是rawdata
                   heatmap_tsne_markers=heatmap_tsne_markers,
                   trans_method="simpleAsinh",
                   groups=groups,
                   #major_cond="Tissue_Type",
                   dot_size=4,
                   #dot_colours = colorRampPalette(c("#4575B4" ,"#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59" ,"#D73027")),
                   show_legend=T,
                   reduction_dm1="UMAP1",
                   reduction_dm2="UMAP2",
                   output_dir="UMAP_heatmap",
                   show_scale=F,
                   show_contour = F,
                   output_format="tiff",
                   single_file=T)

}
      
#===========================13.降维和聚类结果做图（Cytofkit Shinny App）============================================

#  (A) Input：降维结果汇总:
#     1) 参数设置
      dimReducedRes <- list(tsne_result,rumap_result)                 #<-- 降维结果（可以将多次结果汇总）
      names(dimReducedRes)<-c("tsne","UMAP")                    #<-- 降维名称（可以将多次结果汇总，名称与前面结果顺序一致）
      clusterRes<-list(combined_transdata_sampled$PhenoGraph,combined_rawdata_sampled$metacluster) #<-- 聚类结果（可以将多次结果汇总）
      names(clusterRes)<-c("Phenograph","FlowSOM")                 #<---聚类名称（可以将多次结果汇总，名称与前面结果顺序一致）

#     2) Checkpoint: 
      str(dimReducedRes)
      str(clusterRes)


#  (B) 生成Cytofkit_ShinnyApp 识别的RData文件
      export_cytofkit_RData(expressionData =combined_transdata_sampled,
                            dimReducedRes = dimReducedRes,
                            clusterRes = clusterRes, 
                            projectName = "cytofkit", #<-- 设置Rdata的文件名
                            rawFCSdir = raw_fcs_dir,
                            resultDir = raw_fcs_dir,
                            dimRedMarkers = colnames(tSNE_input_data), #<-括号内设为降维分析使用的表达矩阵数据
                            sampleNames = unique(File_ID))

##Checkpoint 如看结果，请新打开一个Rstudio窗口，打开Open_ShinnyApp.R 运行里面的语句。


#===========================14.将结果导出成FCS文件==================================================================================

#  (A) Input： 整合待导出数据

     combined_data_output <- combined_rawdata_sampled[,0]
     if(Do_PhenoGraph)  combined_data_output <- data.frame(combined_data_output,PhenoGraph = combined_rawdata_sampled$PhenoGraph)  
     if(Do_FlowSOM)     combined_data_output <- data.frame(combined_data_output,metacluster = combined_rawdata_sampled$metacluster) 
     if(Do_tSNE)        combined_data_output <- data.frame(combined_data_output,tsne_result)  
     if(DO_UMAP)        combined_data_output <- data.frame(combined_data_output,rumap_result) 

          
#  (B) CheckPoint:

     head(combined_data_output)

#  (C) 根据File_ID将合并数据还原成单个样本数据，导出FCS文件

     cytof_addToFCS_modified(
        data=combined_data_output,
        rawFCSdir=raw_fcs_dir,
        analyzedFCSdir=paste0(wdir,"/tsne_UMAP_output"),
        newHeader = "tSNE_UMAP_"  #<-- 可以自行设定输出文件的前缀
        )
#==========================END============================================================================================
