
################################################################## #
#                                                                  #
#    版本信息：visual_cluster_stat           发布时间：2021-03-30  #
#    Pipeline versiion: 2.10                                       #
#                                                                  #
################################################################## # 




#======================================1.载入需要的工具包==============================================

rm(list = ls())#清除内存变量
gc()
library(cytofkit)
library(dplyr)
library(gplots)
library(ggrepel)
library(Rmisc)
library(colorRamps)
library(RColorBrewer)
library(car)
library(ggpubr)
library(cytofexplorer)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(egg)
library(flowCore)
library(DescTools)

#======================================2. 目录设置==============================================================

#  (A) Input：<--
        raw_fcs_dir="21_Phenograph_tsne/output"  #<---设置读取数据的文件夹名称

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
        
#======================================3.数据预处理（读取->Downsample->合并）=================

#3.1读取->Downsample->合并

#  (A) Input：<--
        mergeMethod <- "all"      #合并方法："ceil", "all", "fixed", "min"
        fixedNum <- 50000         # 设置从每个文件抽取的细胞数目

#  (B) 读取文件，DownSample，Merge
        file_name <- list.files(raw_fcs_dir,pattern='.fcs$', full=TRUE)
        combined_data_raw <- cytof_exprsMerge(fcsFiles = file_name,
                                              transformMethod = "none",
                                              mergeMethod =mergeMethod,
                                              fixedNum = fixedNum)
#  (C) 简化列名
        paraname<-colnames(combined_data_raw)
        paraname<-sub("<NA>","",paraname)
        paraname<-sub(".*<","",paraname)
        paraname<-sub(">.*","",paraname)
        paraname<-sub("-","_",paraname)
        paraname<-sub("^[0-9][0-9][0-9][A-Z][a-z]_|^89Y_","",paraname)  #去掉数字标签
        paraname<-sub("_EQ[0-9]$","",paraname)  #去掉EQ标签
        #paraname<-sub("+","_",paraname)
        colnames(combined_data_raw)<-paraname
#  (D) 增加File_ID
        File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
        combined_data_raw<-data.frame(combined_data_raw,File_ID)


#  (E) Checkpoint
head(combined_data_raw)


#===========================4. 追加MetaData数据========================================================

#  (A)  marker选择
#       打开all_markers.csv，增加expr_para、heatmap两列分别用来指定用来分析差异表达和出现在heatmap中的marker,
#       设置好通道后，save

        all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),header = TRUE,stringsAsFactors = F)
        markernames<-colnames(combined_data_raw)
        markernames[c(1:(nrow(all_markers)-1))]<-all_markers[-nrow(all_markers),"markers"]
        colnames(combined_data_raw)<-markernames

#  (B)  分组信息读取
        groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors=FALSE)
        groups$File_ID<-unique(File_ID)

#  (C)checkpoint: -->  查看groups内容,确认分组信息被正常读取
        groups
        all_markers

        

#======================================5.数据初步整理，全局参数设定=========================================================
 
cluster_stat<- stat_by_cluster( combined_data_raw,
                      all_markers,
                      #cluster_merge = list(c(0)),
                      cluster_name="PhenoGraph",  #<——要进行统计的cluster的通道名称
                      alias_as_label=F,
                      summerise_method="median",   #<--均值统计的方法："median"(中位数)或者"mean"(均值)
                      groups=groups,
                      major_cond="Group_1",    #<--all_samples.csv的列名，指明要寻找差异的分组形式
                      group_seq=c("HC","MMD"), #设置major_cod 在各个统计图中的顺序
                      stat.paired=F,               #<--是否为成对样本
                      stat.method="wilcox-test"         #<--统计的方法："t-test"或者"wilcox-test"
                      )

#======================================6.亚群丰度分析barplot,heatmap,histogram=========================================================================
#绘制boxplot
library(RColorBrewer)
library("ggsci")
abundance_metadata <-draw_abundance_report(
                     cluster_stat,
                     #cluster_id=c(1:10),            #<——指定要输出的cluster id
                     heatmap_ctrl=c("HC"),       #<--指定heatmap中的control组
                     stat.paired=F,                  #<--是否为成对样本

                     #heatmap图形设置：
                     #heatmap1_color =colorRampPalette(c("black","yellow")),
                     #heatmap2_color =colorRampPalette(c("blue","white","red")),
                     Rowv=T,Colv=T,dendrogram="both",#<--heatmap 行列聚类和树形图设置
                     
                     #反映各个样本亚群组成的barplot的设置：
                     barplot_clustered =T,          #<--样本之间是否进行聚类
                     barplot_direction="v",          #<--barplot的方向是横向"h"还是纵向"v"
                     
                     #在boxplot上显示各组两两之间的p值 
                     comparisons = list(c("HC","MMD"))
                     #comparisons.stat.paired=T
                     )

#======================================7.亚群丰度分析volcano plot=========================================================
#(此段代码可以复制，输出不同组之间的火山图)
draw_abundance_volcano(cluster_stat,
                       cond1="HC",       #火山图左侧的组
                       cond2="MMD")     #火山图右侧的组


#======================================8.Marker表达差异分析(barplot,heatmap)===========================================
##需要all_markers.csv中含有expr_para一列
#生成每个cluster的expr_para的heatmap,boxplot
expr_metadata<-cluster_expr_report(
                      cluster_stat,
                      #cluster_id=c(1,8),                #<——指定要输出的cluster id
                      heatmap_ctrl=c("HC"),          #<--指定heatmap中的control组
                      stat.paired=F,                     #<--是否为成对样本
                      #hide_ctrl = T,
                      
                      #heatmap图形设置：
                      #heatmap1_trans_method="0_to_Max",
                      #heatmap1_color =colorRampPalette(c("black","yellow")),
                      #heatmap2_color =colorRampPalette(c("blue","white","red")),
                      Rowv=T,Colv=T,dendrogram="both",   #<--heatmap 行列聚类和树形图设置
                      
                      #在boxplot上显示各组两两之间的p值 
                      comparisons = list(c("HC","MMD")),
                      comparisons.stat.paired=F
                      )

#======================================9.Marker表达差异分析volcano plot========================================================
#绘制差异表达火山图(此段代码可以复制，输出不同组之间的火山图)
draw_expr_volcano(cluster_stat,
                  cond1="HC",           #火山图左侧的组
                  cond2="MMD",         #火山图右侧的组
                  color_by_significant=T,
                  label_repel=F
                  )         

#======================================10.绘制heatmap和density map==============================================
##需要all_markers.csv中含有heatmap一列

draw_expr_heatmaps(combined_data_raw,
                   cluster_stat,
                   #cluster_id=c(1:5),
                   facet_by_groups=F,                 #是否分开显示不同组的heatmap
                   Rowv=T,
                   Colv=T,
                   trans_method="simpleAsinh",        #数据转换函数
                   #max_b=0.5,
                   #是否在heatmap里面显示abundance barplot
                   show_abundance_barplot=T,
                   
                   #denstiy plot的颜色设置,保持默认即可
                   output_density_plot=T,
                   density_plot_fill_by="expression",
                   density_plot_color_by="none",
                   density_bkgd_fill_by="expression"
                   #groups_to_show=c("PBMC","Biopsy") #指定显示哪些数据
                   )


