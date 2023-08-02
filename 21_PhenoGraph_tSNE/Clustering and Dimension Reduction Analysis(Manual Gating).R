
################################################################### #
#                                                                   #
#  版本信息：clustering and DR(manual_gating)  发布时间：2020-3-30  #
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
library(DescTools)

#===========================1.Pipeline设置==========================================================

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
      markers2csv(paraname_raw,colnames(combined_data_raw))

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
      
#===========================6a. 生成Downsample的FCS文件=====================================================

#  (A) Input： 整合待导出数据

      combined_data_output <- combined_data_raw[,0]
      combined_data_output <- data.frame(combined_data_output,downsample = 1)  #<-- 输入待导出的通道


#  (B) CheckPoint:
      head(combined_data_output)

#  (C) 根据File_ID将合并数据还原成单个样本数据，导出FCS文件
      cytof_addToFCS_modified(
        data=combined_data_output,
        rawFCSdir=raw_fcs_dir,
        analyzedFCSdir=paste(wdir,"/downsample",sep = ""),
        newHeader = "" #<-- 可以自行设定输出文件的前缀
      )

#===========================6b. 从XML File中提取cellular gating ref和 gating strategy信息（Cytobank版） ================================================

#本模块的目的是生成每个文件的gating index文件以及所有population的gating_strategy文件\

library(CytoML)
library(flowWorkspace)

#  (A)整合xml文件和fcs文件数据，得到包含单细胞数据和gating数据的gateset object：gs      
     file_name <- list.files(pattern='.fcs$',path = './downsample', full=TRUE)
     xmlfile   <- list.files(pattern='.xml', path=raw_fcs_dir, full.names=TRUE)
     gs <- cytobank2GatingSet(xmlfile, file_name) 

#  (B)从gs中提取单细胞gating数据，输出为每个FCS的cell gating ref文件
     generate_gating_ref(gs)

#  (C)从gs中提取每个population的manual gating strategy文件
     generate_gating_strategy(gs)
                                                                          
#  (D)checkpiont:(运行下面代码观察相关文件是否正常输出)
     cat(list.files('./downsample',pattern='.csv$', full=F))
     cat(list.files(metadata_dir,pattern='.csv$', full=F))

#===========================6c.Extract manual gate population results========================================================

#  (A)读取以及合并gating_ref文件：      
     gating_ref_csv <- list.files('./downsample',pattern='.csv$', full=TRUE)
     combined_gating_ref <- gating_ref_merge(gating_ref_csv)

#  (B)提取单细胞的population数据:
#     1) Input<--
#      到metadata目录中找到gating_strategy.csv, 选择需要提取的population：
#      在Pop_Extract_Seq一列输入1,2,3,4...，程序会从小到大的顺序处理各个population。
#      注意：如果一个细胞同时属于两个以上的population，后面population index会覆盖前面的。

     gating_strategy_table<-read.csv(paste0(metadata_dir,"/gating_strategy.csv"),stringsAsFactors = F,header = T)
     manual_gate_result<-extract_pop_ref(combined_gating_ref,gating_strategy_table)

#  (C)checkpiont 
     head(manual_gate_result)
     

#===========================7.（*可选步骤）:手动合并亚群(Manually merge clusters)==============================================

    if(FALSE){    #<--如果要进行合并括号中改为TRUE
    
      merge_result<-cluster_merge(cluster_result=PhenoGraph_result, #<--待合并的聚类结果  
                    merge_list=list(c(1,2,3),c(6:9)),               #<--指出要合并的亚群例如：  list(c(1,2,3),c(6:9))
                    cluster_rename=T)                               #<--是否对cluster重新命名
      
      #Checkpoint
      hist(merge_result,unique(merge_result))
      }

#===========================8. Cluster-marker表达Preview(Density和Heatmap)====================================================
##整合数据，从各个Meta cluster或文件中抽取定量细胞进行可视化(Downsample)

#  (A) Input: 整合数据：
      combined_rawdata_analysed <- data.frame(combined_data_raw,
                                              manual_gate=manual_gate_result) #<-- 加入待整合的聚类结果
      
#  (B) 生成Marker-Cluster Heatmap和Density图片

      cluster_marker_preview(combined_rawdata_analysed,
                             groups = groups,
                             all_markers=all_markers,
                             cluster_color=dif_seq_rainbow,  #选择颜色 其他选择：rainbow，dif_seq_rainbow,brewer_color_sets
                             cluster_name="manual_gate",      #cluster 通道的名称
                             colored_by="expression"         #"expression"
                             #cluster_id=NULL                # 选择要显示的cluster
                             )
      
      head( combined_rawdata_analysed)
      
#===========================9. 将结果导出成FCS文件=====================================================
      
#  (A) Input： 整合待导出数据

      combined_data_output <- combined_data_raw[,0]
      combined_data_output <- data.frame(combined_data_output,manual_gate = manual_gate_result)  #<-- 输入待导出的通道
     
      
#  (B) CheckPoint:
      head(combined_data_output)
      
#  (C) 根据File_ID将合并数据还原成单个样本数据，导出FCS文件
      cytof_addToFCS_modified(
        data=combined_data_output,
        rawFCSdir=raw_fcs_dir,
        analyzedFCSdir=paste(wdir,"/output",sep = ""),
        newHeader = "manualgate_" #<-- 可以自行设定输出文件的前缀
        )

#===========================10. DownSample==============================================================  

#  (A) Input： 整合rawdata和聚类结果
      combined_rawdata_analysed <- combined_data_raw
      combined_rawdata_analysed <- data.frame(combined_rawdata_analysed,
                                              manual_gate = manual_gate_result)  #<-- 输入待导出的通道
      
      

#  (B) downsample         
      combined_rawdata_sampled<-equal_sample(x=combined_rawdata_analysed,
                                          sample_method="ceil",         #<-- 两种取值："ceil", "all"
                                          sample_type="by_file",        #<-- "by_cluster","by_file","by_cond"
                                          sample_num=200,              #<-- 每个cluster/file/cond抽取的细胞数
                                          cluster_name="manual_gate",    #<-- Cluster的名称
                                          groups=groups,                #仅在by_cond时需要，其他时候可以注释掉
                                          sample_cond = "Tissue_Type"   #仅在by_cond时需要，其他时候可以注释掉
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
#===========================11b. 降维结果可视化========================================================
if(Do_tSNE){

#  (A) Input：
#    1) 将rawdata和降维结果整合：
      combined_rawdata_plot   <- data.frame(combined_rawdata_sampled,
                                      tsne_result)           #<-- 降维结果
#  (B) Checkpoint:
      head(combined_rawdata_plot)
      #head(combined_transdata_plot)

#  (C) 生成tsne-Phenograph系列图片
      draw_tsne_figs(combined_rawdata_plot,  
                     groups=groups,
                     cluster_color=dif_seq_rainbow,
                     cluster_name="manual_gate",
                     alias_as_label = T,
                     major_cond="Tissue_Type",
                     cluster_labels_size=15,
                     dot_size = 3,
                     show_cluster_labels = T,
                     reduction_dm1="tsne_1",
                     reduction_dm2="tsne_2",
                     show_contour=T,
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
                         single_file=T)

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
               cluster_name="metacluster",
               major_cond="Tissue_Type",
               cluster_barplot_data=metacluster_result,
               cluster_labels_size=15,
               dot_size = 3,
               show_cluster_labels = T,
               reduction_dm1="UMAP1",
               reduction_dm2="UMAP2",
               output_dir="UMAP",
               show_contour=T,
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
                   major_cond="Tissue_Type",
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


#===========================14.将结果导出成FCS文件================================================================

#  (A) Input： 整合待导出数据

     combined_data_output <- combined_rawdata_sampled[,0]
     if(DO_UMAP)        combined_data_output <- data.frame(combined_data_output,rumap_result) 
                        combined_data_output <- data.frame(combined_data_output,pop_index = combined_rawdata_sampled$pop_index) 
      
          
#  (B) CheckPoint:

     head(combined_data_output)

#  (C) 根据File_ID将合并数据还原成单个样本数据，导出FCS文件

     cytof_addToFCS_modified(
        data=combined_data_output,
        rawFCSdir=raw_fcs_dir,
        analyzedFCSdir=paste0(wdir,"/tsne_UMAP_output"),
        newHeader = "tSNE_UMAP_"  #<-- 可以自行设定输出文件的前缀
        )

#====================================END===============================================================
