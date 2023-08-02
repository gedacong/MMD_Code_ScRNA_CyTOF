
#2019_10_28 combined_stat函数引入错误捕获机制(try)，计算出错时直接输出warning信息，不会中止
#2019_11_4  在cluster_expr_report里面加入groups_to_show 参数,引入错误捕获机制，防止出图错误中断程序
#2019_11_26 修复cluster_expr_report中的 color_cond bug
#2019_12_5  修复单个subgroup时，没有color_cond的错误
#2019_12_19 在expr_report中增加Rowv=T,Colv=F,dendrogram="row",几个heatmap 聚类参数
#2019_12_20 对cluster_merge的功能进行完善,支持将所有cluster合并为一个进行分析(即对整个文件进行统计)
#2019_12_20 对clusterlabel进行改进,保持合并后cluster与label的对应关系
#2020_01_22 对expr-report进行改版，提高代码效率，取消去除极值处理，增加hide_ctrl参数
#2020_01_30 对expr_report hide_ctrl参数进行扩展（可以同时导出带有ctrl和不带有ctrl的figure）
#2020_09_04 matches-->one_of
#2020-09-04 修复na.omity引起的bug：FDR_pValue <-FDR_pValue[complete.cases(FDR_pValue[,c("p_values","p_adjust","log_ratio","log_p_values")]),]
#2020-10-3  cluster_expr_report函数中增加heatmap1_trans_methodc参数

#' @title combined_stat 统计分析函数
#'
#' @description
#' 对进行单组或两组数据进行统计分析，并对结果进行adjust。也可以为数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'
#' @param cond1_data      一个包含数字的向量（vector），待统计分析数据
#' @param cond2_data      一个包含数字的向量（vector），待统计分析数据，如进行单组分析，cond2_data=NULL
#' @param alternative     一个字符串，反应假设检验的种类,有三种取值："two.sided" (默认), "greater" or "less"；
#' @param mu              一个数字，在单样本本统计时表示平均值的真实值，或者如果您正在进行两个样本测试，则表示均值差异（此时一般取值为0）
#' @param tranform.method 决定计算p值前进行何种transform，“raw” 不做任何变换，"asinh"按照以下公式进行arcsinh转换:y=arcsinh(x/5)；
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @param  silent         逻辑变量，是否输出中间结果
#' @return 返回一个数据框（data.frame),包含p_value,p_adjust,统计方法(stat.method),方差齐性（set.equal.var）信息
#' @export

combined_stat<-function(cond1_data,
                         cond2_data,
                         alternative="two.sided",
                         mu=0,
                         transform.method="raw",
                         stat.paired=FALSE,
                         conf.level=0.95,
                         stat.method="t-test",
                         set.equal.var="FALSE",
                         p.adjust.method="BH",
                         silent=F){
s_methods<-NULL  
s_method<-NULL
p_values<-NULL
var.equals<-NULL
if(stat.method=="t-test"){
  s_method<-"t-test"}else if(stat.method=="wilcox-test"){
    s_method<-"wilcox-test"}

if( transform.method=="asinh" ){ 
  cond1_data<-apply(cond1_data, 2,simpleAsinh)
  cond2_data<-apply(cond2_data, 2,simpleAsinh)
  
}



for(i in 1:ncol(cond1_data)){
  if(silent==F) cat(paste0("Treating with cluster",i,":  "))
  var.equal<-NA

  if(stat.method=="auto"){
    #判断是否是正态分布

    shapiro_result1=shapiro.test(cond1_data[,i,drop=T])$p.value
    if(!is.null(cond2_data)){
    shapiro_result2=shapiro.test(cond2_data[,i,drop=T])$p.value}else{
    shapiro_result2=1
    }

    if(shapiro_result1>0.05 & shapiro_result2>0.05){
      s_method<-"t-test" }else{
        s_method<-"wilcox-test"}
  }
  #t-test
  if(s_method=="t-test"){


    #方差齐性检测
    if(set.equal.var=="auto")
    {var.equal=FALSE
    if(!is.null(cond2_data)){

    pool_data<-c(cond1_data[,i,drop=T],cond2_data[,i,drop=T])
    pool_group<-as.factor(c(rep(1,length(cond1_data[,i,drop=T])),rep(2,length(cond2_data[,i,drop=T]))))
    equal_var<-leveneTest(y =pool_data, group=pool_group)
    if(equal_var$`F value`[1]<0.05) {
      var.equal=FALSE }else{
        var.equal=TRUE
      }
    }
    }else if(set.equal.var=="TRUE"){
             var.equal=TRUE} else if((set.equal.var=="FALSE") ) {
               var.equal=FALSE}

    cond1_with_na<- length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T])
    cond2_with_na<- suppressWarnings(length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T]))
    if(cond1_with_na & cond2_with_na ){
      pv<-try(t.test(x=cond1_data[,i,drop=T],
                 y=cond2_data[,i,drop=T],
                 alternative = alternative,
                 paired = stat.paired,
                 var.equal=var.equal,
                 mu=mu,
                 conf.level=conf.level),
                 silent = T#$p.value
            )
      if ('try-error' %in% class(pv)) 
      {pv<-NA
      message(paste0("warning: Failed to calculate p value of cluster",i,"\n"))
      cat("cond1 data is: ",cond1_data[,i,drop=T],"\n")
      cat("cond2 data is: ",cond2_data[,i,drop=T],"\n")
      } else{pv<-pv$p.value}
    }else{
      pv<-NA
    }
    #p_values<-cbind(p_values,pv)

  }
  else if(s_method=="wilcox-test"){
    #if(nrow(na.omit(cond1_data))>1 &  nrow(na.omit(cond2_data))>1){
    #if(length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T]) &  length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T])){
    cond1_with_na<- length(na.omit(cond1_data[,i,drop=T]))==length(cond1_data[,i,drop=T])
    cond2_with_na<- suppressWarnings(length(na.omit(cond2_data[,i,drop=T]))==length(cond2_data[,i,drop=T]))
    if(cond1_with_na & cond2_with_na ){
    pv<-try(wilcox.test(cond1_data[,i,drop=T],
                      cond2_data[,i,drop=T],
                      alternative = alternative,
                      stat.paired = stat.paired,
                      paired = stat.paired,
                      mu=mu,
                      conf.level=conf.level),
                      silent=T#$p.value
           )
    if ('try-error' %in% class(pv)) 
          {pv<-NA
           message(paste0("warning: Failed to calculate p value of cluster",i,"\n"))
           cat("cond1 data is: ",cond1_data[,i,drop=T],"\n")
           cat("cond2 data is: ",cond2_data[,i,drop=T],"\n")
           } else{pv<-pv$p.value}
    
    }else{

      pv<-NA
    }
    #p_values<-cbind(p_values,pv)

  }

  if(silent==F) cat(paste0(" ",s_method," var.equal=",var.equal,"\n"))
  #p_values<-cbind(p_values,pv)
  #s_methods<-rbind(s_methods,s_method)
  #var.equals<-rbind(var.equals,var.equal)

  p_values<-c(p_values,pv)
  s_methods<-c(s_methods,s_method)
  var.equals<-c(var.equals,var.equal)


}

#adjust p value
p_adjust<-p.adjust(p_values,p.adjust.method,length(p_values))

stat_result<-data.frame(p_values=p_values,
                        p_adjust=p_adjust,
                        s_methods=s_methods,
                        var.equals=var.equals)

return(stat_result)

}



#' @title Statistics of clusters(Abundance and Experssion)
#'
#' @description
#' 统计不同分组中各个cluster中marker表达的差异
#'
#' @param combined_data_raw  在Pepline前面部分读入并合并的表达数据
#' @param all_markers        从 “all_markers.csv”读入的dataframe，包含marker的选择信息
#' @param cluster_name       统计所依据的聚类参数名称，例如"PhenoGraph"、"metacluster"
#' @param alias_as_label     逻辑变量，决定是否使用cluster的别名做为其label


#' @param cluster_merge      指定要合并的组别例如：list(c(1:3),c(5,8)) 意思是将1~3合并，5和8合并

#' @param summerise_method   统计过程中使用的合计方法：有"mean"或者 "median"两种
#' @param groups             groups 从“all_sample.csv”中读入的dataframe，包含实验分组信息（也就是元数据）
#'
#' @param major_cond         在groups所提供的各种分组中，选择一种做为统计差异的依据，例如"Tissue_Type"
#' @param group_seq          设置各组顺序，格式实例： group_seq=c("PBMC","Biopsy")，注：括号内为各组名称；

#'下面5个参数是为后面统计分析中设置的默认值

#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"

#' @return a List containing the data which needed in drawing of heatmap, boxplot and volcano chart
#' @export

stat_by_cluster<- function (combined_data_raw,
                            all_markers,
                            cluster_name="metacluster",
                            cluster_merge=NULL,
                            alias_as_label=FALSE,
                            summerise_method="median", #"mean" or "median" the static method used in summerise
                            groups,

                            #group 参数
                            major_cond,
                            group_seq=NULL,

                            #全局统计默认参数设置
                            stat.paired=FALSE,
                            stat.method="t-test",
                            set.equal.var="FALSE",
                            conf.level=0.95,
                            p.adjust.method="BH"

                          ){

 if(0){
   
   cluster_name="manual_gate"
   cluster_merge=NULL
   summerise_method="median" #"mean" or "median" the static method used in summerise
   groups
   alias_as_label=TRUE
   
   #group 参数
   major_cond
   group_seq=NULL
   
   
   #全局统计默认参数设置
   stat.paired=FALSE
   stat.method="t-test"
   set.equal.var="FALSE"
   conf.level=0.95
   p.adjust.method="BH"
}
  # cluster_name="metacluster"
  # summerise_method="median" #"mean" or "median" the static method used in summerise



cluster_name<-as.character(cluster_name)
ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)
trans_markers<-as.character(subset(all_markers,transform==1)$markers)
used_markers<-union(ht_markers,trans_markers)

cond_name<-major_cond
cond_name<-as.character(cond_name)


if (!is.null(group_seq)){
  groups[,cond_name]<-factor(groups[,cond_name,drop=T],
                             levels=group_seq,
                             ordered=T)
}



if(alias_as_label==TRUE){
  cluster_label_df<-get_cluster_alias(cluster_name,raw_fcs_dir)
  #cluster_levels<- levels(cluster_label_df[,2])
  #cluster_label_df[,2]<-as.character(cluster_label_df[,2])
  combined_data_raw<-left_join(combined_data_raw,cluster_label_df,by=cluster_name)
  cluster_name<-colnames(cluster_label_df)[2]
  }
  
combined_data_df<- tbl_df(combined_data_raw)

if(!is.null(cluster_merge)){
  for(merge_cluster_i in c(1:length(cluster_merge))){
    #merge_cluster_i=1
    merge_clusters<-paste(cluster_merge[[merge_cluster_i]],collapse=", ")
    merge_cluster_name<-max(combined_data_df[,cluster_name,drop=TRUE])+1
    message(paste0("Note --Cluster",merge_clusters," has been merged as cluster ",merge_cluster_name,"\n",collapse=","))
    merge_cluster_id<-which(combined_data_df[,cluster_name,drop=TRUE] %in% cluster_merge[[merge_cluster_i]])
    combined_data_df[merge_cluster_id,cluster_name]<-merge_cluster_name
   }
}

  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }

  #generating data for total heatmap

  cluster_expr_matrix<-
    combined_data_df %>%
    group_by_at(cluster_name) %>%
    summarise_if(is.numeric,use_method,na.rm=TRUE)
  

  combined_data_tidy<-
    combined_data_df %>%
    group_by_at(c("File_ID",cluster_name)) %>%
    summarise_if(is.numeric,use_method,na.rm=TRUE)

  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)
  cluster_id<-which(colnames(combined_data_tidy)==cluster_name)
 
  check<- function(var) {
          return(var %in%  unique(combined_data_tidy[,cluster_id,drop=T]))
    
      }
  #Pop_names<-Manual_Pop_Index[as.vector(unlist(lapply(Manual_Pop_Index[,"Pop_index",drop=T], check))),"Pop_name"]
  
  clustern<-nrow(unique(combined_data_tidy[cluster_id]))

  cluster_expr_trans<-list()
  cluster_expr_ctrl<-list()
  cluster_stat<-list()
  cluster_cell_count<-list()
  cluster_expr_raw<-list()

  
  
  merged_cluster_id<-sort(unique(combined_data_tidy[,cluster_id,drop=T]))
  
  for(cluster_i in merged_cluster_id){
    
  # cluster_i=7
    #cluster_i=unique(combined_data_tidy[,cluster_id,drop=T])[1]
    cat(paste0("Start to treat with cluster: ",cluster_i,"\n"))

    #统计指定cluster上的表达数据
    combined_data_tidy$File_ID<-as.character(combined_data_tidy$File_ID)
    # expr_raw<-
    # combined_data_tidy %>%
    #       dplyr::filter_(paste0(cluster_name,"==",cluster_i)) %>%
    #       full_join(groups,by="File_ID") %>%
    #       arrange_("Short_name") %>%
    #       arrange_(cond_name)
   
    expr_raw<-
    combined_data_tidy %>%
      dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.==cluster_i)) %>%
      full_join(groups,by="File_ID") %>%
      arrange_at("Short_name") %>%
      arrange_at(cond_name) 
    expr_raw$metacluster<-cluster_i #去除NA
    
#str(cluster_name)
    #which(is.na(combined_data_tidy$metacluster))
    #which(is.na(expr_raw$metacluster))
    #统计指定cluster中含有各个文件的细胞数目
    # 
    combined_data_df$File_ID<-as.character(combined_data_df$File_ID)
    #cell_count<-
    # combined_data_df %>%
    #   dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.==cluster_i)) %>%
    #       count_("File_ID") %>%
    #       full_join(groups,by="File_ID") %>%
    #       arrange_at("Short_name") %>%
    #       arrange_at(cond_name)
    # 
    # 
    cell_count<-  
    combined_data_df %>%
      dplyr::filter_at(vars(one_of(cluster_name)),all_vars(.==cluster_i)) %>%
      group_by_at("File_ID")%>%
      tally() %>%
      full_join(groups,by="File_ID") %>%
      arrange_at("Short_name") %>%
      arrange_at(cond_name)
    
    
    
    expr_raw<-ungroup(expr_raw)


    #计算Arcsin ratio
    expr_trans<-expr_raw
    expr_trans[,trans_markers] <- apply(expr_trans[, trans_markers,drop = FALSE],2,simpleAsinh)
    expr_ctrl<-as.data.frame(matrix(0,ncol(expr_trans),nrow=0))
    colnames(expr_ctrl)<-colnames(expr_trans)
    #which(is.na(expr_trans$metacluster))
    

    #col_color_index<-as.numeric(as.factor(t(expr_raw[,cond_name,drop=T])))
    cluster_expr_trans[[cluster_i]]<-expr_trans
    cluster_expr_ctrl[[cluster_i]]<-expr_ctrl
    cluster_cell_count[[cluster_i]]<-cell_count
    cluster_expr_raw[[cluster_i]]<-expr_raw

    
    #which(is.na(cluster_expr_trans[[7]]$metacluster))
    
  }


  #Abundance数据汇总
  cat("Summarising cluster abundance data: Cell Count and Percentage\n ")
  cluster_count_summary<-cluster_cell_count[[merged_cluster_id[1]]][,1]
  for(cluster_i in merged_cluster_id)
   # cluster_i<-merged_cluster_id[1]
  { cluster_count_summary=full_join(cluster_count_summary,cluster_cell_count[[cluster_i]][,c(1,2)],by="File_ID")
  }

  cluster_count_summary[is.na(cluster_count_summary)]<-0
  
  
  
  # if(is.null(Manual_Pop_Index)){
  
   if (alias_as_label==F){
            colnames(cluster_count_summary)<-c("File_ID",paste0("cluster",merged_cluster_id))
    } else{
            colnames(cluster_count_summary)<-c("File_ID",as.character(merged_cluster_id) )   
             }
  # }else{
  #           colnames(cluster_count_summary)<-c("File_ID",Manual_Pop_Index$Pop_name)
  #           names(cluster_count_summary)<-c("File_ID",Manual_Pop_Index$Pop_name)
  #           }
  # 
  percentage<-function(var){return(var/sum(var)*100)}
  
  cluster_percent_summary<-data.frame(cluster_count_summary,check.names = F)
  
  
  if (ncol(cluster_percent_summary)==2) {
  cluster_percent_summary[,-1]<-as.matrix(apply(cluster_percent_summary[,-1,drop=F],1,percentage))} else{
    
  cluster_percent_summary[,-1]<-t(apply(cluster_percent_summary[,-1,drop=F],1,percentage))  
  }
  
  
  
  
  cluster_count_summary$File_ID<-as.character(cluster_count_summary$File_ID)
  cluster_percent_summary$File_ID<-as.character(cluster_percent_summary$File_ID)



  cluster_stat[["metadata"]]<- data.frame( cluster_name=cluster_name,
                                           clustern=clustern,
                                           alias_as_label=alias_as_label,
                                           summerise_method=summerise_method,
                                           #heatmap 参数
                                           major_cond=major_cond,
                                           #heatmap_ctrl=heatmap_ctrl,
                                           stat.paired=stat.paired,
                                           stat.method=stat.method,
                                           p.adjust.method=p.adjust.method,
                                           set.equal.var=set.equal.var,
                                           conf.level=conf.level
                                           )
  # if(!is.null(Manual_Pop_Index)){
  # Manual_Pop_Index$cluster<-paste0("cluster",Manual_Pop_Index$Pop_index)}
  
  
  cluster_stat[["cluster_expr_matrix"]]<-cluster_expr_matrix
  cluster_stat[["cluster_expr_trans"]]<-cluster_expr_trans
  #cluster_stat[["cluster_expr_ctrl"]]<-cluster_expr_ctrl
  #cluster_stat[["volcano_data_all"]]<-volcano_data_all
  #cluster_stat[["col_color_index"]]<-col_color_index
  cluster_stat[["cluster_cell_count"]]<-cluster_cell_count
  cluster_stat[["ht_markers"]]<-ht_markers
  cluster_stat[["groups"]]<-groups
  cluster_stat[["all_markers"]]<-all_markers
  cluster_stat[["cluster_count_summary"]]<-cluster_count_summary
  cluster_stat[["cluster_percent_summary"]]<-cluster_percent_summary
  cluster_stat[["cluster_expr_raw"]]<-cluster_expr_raw
  cluster_stat[["group_seq"]]<-group_seq
  cluster_stat[["merged_cluster_id"]]<-merged_cluster_id
  if (alias_as_label==T){
      cluster_stat[["cluster_label_df"]]<-cluster_label_df
  }
  # which(is.na(cluster_stat[["cluster_expr_trans"]][[7]]$metacluster))
  # which(is.na(cluster_expr_trans[[7]]$metacluster))
  return(cluster_stat)
}






#'@title 生成各个cluster的Heatmap以及boxplot
#'@description
#'依据stat_by_cluster生成的list数据，生成
#'1）heatmap1：指定cluster中，各个样本中的marker的表达情况
#'2）heatmap2：指定cluster中，各个样本中的marker与对照的相对强度（Arcsinh raitio）
#'3）boxplots：统计各marker在该cluster中的表达
#'
#'@param cluster_stat  stat_by_cluster()函数生成的list数据
#'@param cluster_id    指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param color_cond     指定boxplot依据哪个条件进行着色
#'@param subgroup_cond      指定heatmap_ctrl对应的分组 
#'@param group_colorset boxplot采用的配色方案：默认 brewer_color_sets, 也可以选择其他公开生成色阶的函数例如 rainbow等  
#'@param output_dir    输出数据文件夹名称

#'@param heatmap_ctrl   在major_cond中选择一个（或多个）在Heatmap中显示做为control的组名，例如"A_Tumor"
#'@param heatmap1_color 设置Arcsinh expression heatmap的colorbar
#'@param heatmap1_trans_method 数据转化方式，有五种：
#'                            "CytofAsinh"，CytofAsinh转换所有表达数据；
#'                            "simpleAsinh",simpleAsinh转换所有表达数据,（公式为asinh(x/5)）,这是默认设置；
#'                            "0_to_Max"，所有Marker的信号强度除以最大值，线性转换，通过除以各通道信号的最大值把数值scale到0~1；
#'                            "Min_to_Max"，线性转换，最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
#'                            "simpleAsinh_0_to_Max" 先Arcsinh转换，然后除以各通道信号的最大值把数值scale到0~1范围内；
#'                            "simpleAsinh_Min_to_Max" 先Arcsinh转换，然后最小值转换成0，最大值转换成1，最大限度展现population的表达差异；
#'
#'
#'@param heatmap2_color 设置Arcsinh ratio heatmap的colorbar
#'@param Rowv,Colv         逻辑变量，分别heatmap设定行和列是否聚类
#'@param dendrogram        显示heatmap行或者列的树形图,"both","row","column","none"
#'@param show_pvalues  逻辑变量，TRUE或者FALSE，决定是否在图上显示p Value
#'@param hide_ctrl         在图片中隐藏对照组（仅用于多个ctrl的情况）


#'@param comparisons    list变量，用来指定需要显示p值的组别，例如list(c("PBMC","Biopsy")，c("PBMC","Tumor"))就是要分别显示PBMC与Biopsy和Tumor两组之间的p值；如果comparisons=NULL，则显示各组与对照组的p值
#'@param comparisons.stat.paired     逻辑变量，TRUE或者FALSE，表示comparisons是否是成对检验
#'@param comparisons.stat.method  指定组别的p值计算方法


#' 
#'如使用统计默认参数，以下参数无需设置
#'@param groups_to_show 想要展示部分组时，指定展示组的组名，默认为NULL，展示所有组
#'@param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#'@param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#'@param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#'@param conf.level      显著性水平，默认0.95,（即p<0.05为显著）
#'@param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#'@return none
#'@export


cluster_expr_report<-function(cluster_stat,
                              cluster_id=NULL,
                              color_cond=NULL,
                              subgroup_cond=NULL,
                              group_colorset =brewer_color_sets,
                              output_dir="Cluster_expr_report(heatmap and boxplot)",
                              
                              #heatmap parameters
                              heatmap_ctrl=NULL,
                              heatmap1_trans_method="simpleAsinh",
                              heatmap1_color =colorRampPalette(c("black","yellow")),
                              heatmap2_color =colorRampPalette(c("blue","white","red")),
                              Rowv=T,Colv=F,dendrogram="row", #heatmap 聚类参数
                              show_pvalues=TRUE,
                              hide_ctrl=NULL,

                              #boxplot parameters
                              comparisons=NULL,
                              comparisons.stat.paired=NULL,
                              comparisons.stat.method=NULL,
                              
                              
                              #如使用全局统计默认参数，以下参数无需设置
                              groups_to_show=NULL,
                              stat.paired=NULL,
                              stat.method=NULL,
                              set.equal.var=NULL,
                              conf.level=NULL,
                              p.adjust.method=NULL){

  
  if(0){

    cluster_id=NULL
    color_cond=NULL
    subgroup_cond=NULL
    group_colorset =brewer_color_sets
    output_dir="Cluster_expr_report(heatmap and boxplot)"
    
    #heatmap parameters
    heatmap_ctrl=NULL
    heatmap1_trans_method="simpleAsinh"
    heatmap1_color =colorRampPalette(c("black","yellow"))
    heatmap2_color =colorRampPalette(c("blue","white","red"))
    Rowv=T
    Colv=F
    dendrogram="row" #heatmap 聚类参数
    show_pvalues=TRUE
    hide_ctrl=F
    
    #boxplot parameters
    comparisons=NULL
    comparisons.stat.paired=NULL
    comparisons.stat.method=NULL
    
    
    #如使用全局统计默认参数，以下参数无需设置
    #group_seq=NULL,
    groups_to_show=NULL
    stat.paired=NULL
    stat.method=NULL
    set.equal.var=NULL
    conf.level=NULL
    p.adjust.method=NULL
    
  }
  #cat(paste0("./",output_dir))
  
  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  copy_para_csv(output_dir,c("expr_report_boxplot.csv","expr_report_countbar.csv","expr_report_heatmap.csv","expr_report_instruction.pptx"))
  
  paralist_h<-get_para(output_dir,"expr_report_heatmap.csv")
  
  eval(parse(text =paralist_h$items))
  paralist_bar<-get_para(output_dir,"expr_report_countbar.csv")
  eval(parse(text =paralist_bar$items))
  paralist_boxplot<-get_para(output_dir,"expr_report_boxplot.csv")
  eval(parse(text =paralist_boxplot$items))
  
  
  
  col_color_index  <-  cluster_stat[["col_color_index"]]
  #col_color<-brewer.pal(length(unique(col_color_index)),name="Set1")[col_color_index]
  #col_color<-brewer_color_sets(length(unique(col_color_index)))[col_color_index]
  cluster_expr_trans<- cluster_stat[["cluster_expr_trans"]]
  cluster_expr_raw<- cluster_stat[["cluster_expr_raw"]]
  #cluster_expr_ctrl <- cluster_stat[["cluster_expr_ctrl"]]
  cluster_cell_count<- cluster_stat[["cluster_cell_count"]]
  ht_markers        <- cluster_stat[["ht_markers"]]
  summerise_method  <- cluster_stat[["metadata"]]$summerise_method
  cluster_name      <- cluster_stat[["metadata"]]$cluster_name
  clustern          <- cluster_stat[["metadata"]]$clustern
  merged_cluster_id <-cluster_stat[["merged_cluster_id"]]
  major_cond       <-as.character(cluster_stat[["metadata"]]$major_cond)
  groups            <-cluster_stat[["groups"]]
  cond_name<-major_cond
  cond_name<-as.character(cond_name)
  # Pop_names<-cluster_stat[["Pop_names"]]
  # Manual_Pop_Index<-cluster_stat[["Manual_Pop_Index"]]
  
  
#which(is.na(cluster_expr_trans[[7]]$metacluster))

  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }
  
  
  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)
  
  if(is.null(color_cond)) color_cond    <- as.character(major_cond)
  

  if(is.null(groups_to_show)){
    groups_to_show=unique(groups[,major_cond,drop=T])  
  } 
  
  if (!all(heatmap_ctrl %in% groups_to_show))
  { message("Error- heatmap_ctrl should be included in groups_to_show\n")
    return(NULL)
  }  
  
  #确定trans_method
  if(heatmap1_trans_method=="CytofAsinh") trans_fun=cytofAsinh
  if(heatmap1_trans_method=="simpleAsinh") trans_fun=simpleAsinh
  if(heatmap1_trans_method=="0_to_Max") trans_fun=normalize_Zero_Max
  if(heatmap1_trans_method=="Min_to_Max") trans_fun=normalize_min_max
  if(heatmap1_trans_method=="simpleAsinh_0_to_Max") trans_fun=simpleAsinh_Zero_Max
  if(heatmap1_trans_method=="simpleAsinh_Min_to_Max") trans_fun=simpleAsinh_Min_to_Max
  
  
  ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)
  if(length(ht_markers)==0){
    stop("Could not get markers list for diffrential expression analysis, Please check the settings of following column in all_markers.csv: expre_para " )
    }
  trans_markers<-as.character(subset(all_markers,transform==1)$markers)
  used_markers<-union(ht_markers,trans_markers)
  #heatmap_ctrl    <- as.character(cluster_stat[["metadata"]]$heatmap_ctrl)
  if(is.null(cluster_id)){
    cluster_id=merged_cluster_id
  }


  if(stat.method=="auto") {
     stat.method="t.test"
     cat("Warning: Outputing fiures in auto mode, set stat.method to t-test\n" )

  }else { stat.method=sub("-",".",stat.method)}

  if(is.null(comparisons.stat.paired)) comparisons.stat.paired<-stat.paired
  if(is.null(comparisons.stat.method)) {comparisons.stat.method<-stat.method}else{
    (comparisons.stat.method=sub("-",".",comparisons.stat.method))}

  


  
  if (!dir.exists(paste0("./",output_dir,"/",cluster_name,"_csvs"))) {
    dir.create(paste0("./",output_dir,"/",cluster_name,"_csvs"))
  }

  wdir=getwd()

  cat(paste0("\nOutput to folder: ",wdir,"/",output_dir,"\n"))



  major_cond<-as.character(major_cond)
  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)

  #cluster_id<-which(colnames(combined_data_raw)==cluster_name)
  #clustern<-max(unique(combined_data_raw[cluster_id]))
 
   expr_raw<-cluster_expr_raw[[cluster_id[1]]] %>%
    dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))
  #生成合并数据
  expr_ctrl_merge<-expr_raw[0,]
  expr_trans_merge<-expr_raw[0,]
  
  
  for(cluster_i in cluster_id){
    
    # if(!is.null(Manual_Pop_Index)) { Manual_name<-Pop_names[cluster_i]}else{
    #   Manual_name<-""}

    #cluster_i=7
    cat(paste0("Output  heatmaps,barplots,boxplots of cluster ",cluster_i,"...\n"))
    

    
    expr_trans<-cluster_expr_trans[[cluster_i]] %>%
                dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))
    
    expr_raw<-cluster_expr_raw[[cluster_i]] %>%
      dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))
    
    expr_trans_ct<-expr_raw  #customerised tranform heatmap
    expr_trans_ct[,ht_markers]<-apply(expr_raw[,ht_markers],2,trans_fun)
    
    
    
    col_color_index<-as.numeric(as.factor(t(expr_trans[,color_cond,drop=T])))
    col_color<-group_colorset(length(unique(col_color_index)))[col_color_index]
    
    #画heatmap
    #Population_Symbol<-rawToChar(as.raw(cluster_i+64))

    # expr_raw<-cluster_expr_raw[[cluster_i]] %>%
    #   dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))
    # 

    
    expr_ctrl=NULL
    

    if (is.null(subgroup_cond))  {
      #expr_raw$whole<-"Yes"
      expr_trans$whole<-"Yes"
      subgroup_cond<-"whole"
    } else if(subgroup_cond=="whole"){
      #expr_raw$whole<-"Yes"
      expr_trans$whole<-"Yes"
      
    }
 
 
    #计算control
    
    if(length(heatmap_ctrl)!=nrow(unique(expr_trans[,subgroup_cond]))) message("Error:Number of Heatmap_ctrl is not equal to number of subgroup_cond")

    for(subgroup_cond_i in t(unique(expr_trans[,subgroup_cond]))){

    #subgroup_cond_i<-"Term"
    expr_subgroup_cond<-as.data.frame(matrix(0,ncol(expr_trans),nrow=0))
    colnames(expr_subgroup_cond)<-colnames(expr_trans)
  
    expr_trans_cond<-dplyr::filter_(expr_trans,paste0(subgroup_cond,"==subgroup_cond_i"))
    #expr_trans_cond<-dplyr::filter_at(expr_trans,vars(one_of(subgroup_cond)),all_vars(.==subgroup_cond_i))

      
    heatmap_ctrl_i<-intersect(unique(expr_trans_cond[,cond_name,drop=T]),heatmap_ctrl)
      
    #计算paired samples
    if(stat.paired)
          {
          for(conditions in t(unique(expr_trans_cond[,cond_name,drop=T]))){  #
              # expr_ctrl_add<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==conditions"))
              
              expr_ctrl_add<-dplyr::filter_at(expr_trans_cond,vars(one_of(cond_name)),all_vars(.==conditions))
              
              #single_subgroup_cond<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==heatmap_ctrl_i"))
              single_subgroup_cond<-dplyr::filter_at(expr_trans_cond,vars(one_of(cond_name)),all_vars(.==heatmap_ctrl_i))
              
              not_groups_para_id<-which(!(colnames(expr_ctrl_add)%in% colnames(groups)))
              expr_ctrl_add[,not_groups_para_id]<-single_subgroup_cond[,not_groups_para_id]
              expr_subgroup_cond<-rbind.data.frame(expr_subgroup_cond,expr_ctrl_add)}
           } else{
            #unpaired_control_cond<-dplyr::filter_(expr_trans_cond,paste0(cond_name,"==heatmap_ctrl_i"))
            unpaired_control_cond<-dplyr::filter_at(expr_trans_cond,vars(one_of(cond_name)),all_vars(.==heatmap_ctrl_i))
            
            single_subgroup_cond<-unpaired_control_cond[1,]
            single_subgroup_cond[,used_markers]<-summarise_if(unpaired_control_cond[,used_markers],is.numeric,funs(use_method))
            expr_subgroup_cond<-expr_trans_cond
            for(row_i in c(1:nrow(expr_subgroup_cond))){
                expr_subgroup_cond[row_i,used_markers]<-single_subgroup_cond[,used_markers]
                }
          }
    
    if(is.null(expr_ctrl))expr_ctrl<-expr_subgroup_cond
    else{ expr_ctrl<-rbind.data.frame(expr_ctrl,expr_subgroup_cond) }

    }
    
    
    #Generate merge data, used as return value
    expr_ctrl_merge <-rbind.data.frame(expr_ctrl_merge, expr_ctrl)
    expr_trans_merge<-rbind.data.frame(expr_trans_merge,expr_trans)

    expr_ctrl<-
      expr_ctrl %>%
          arrange_at("Short_name") %>%
               arrange_at(cond_name)
    
    expr_trans_ratio<-expr_trans
    expr_trans_ratio[,ht_markers]<-expr_trans[,ht_markers]-expr_ctrl[,ht_markers]
    rownames(expr_trans_ratio)<-expr_trans_ratio$Short_name
    
    expr_report<-function(hide_ctrl){
            file_tag=""    
         if(hide_ctrl==T){
              #message(paste0("hide_ctrl is set to TRUE,heatmap_ctrl groups were hided in all output figures\n"))
              expr_trans<-dplyr::filter_at(expr_trans,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
              expr_trans_ct<-dplyr::filter_at(expr_trans_ct,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
              expr_trans_ratio<-dplyr::filter_at(expr_trans_ratio,vars(one_of(major_cond)),all_vars(!(. %in% heatmap_ctrl)))
              file_tag<-"(HideCtrls)"
              
            }
            
        #======================生成heatmap================================================
            #定义heatmap函数

            
            draw_heatmap<-function(xdata,
                                   heatmap_color,
                                   key.title){
              #xdata<-expr_trans
              #cat("Output heatmap\n")
              expr_heatmap_data<-xdata %>%
                arrange_at("Short_name") %>%
                arrange_at(major_cond)%>%
                data.frame()
              color_cond_names<-expr_heatmap_data[,color_cond,drop=T]
              col_color_index<-as.numeric(factor(color_cond_names,levels = unique(color_cond_names),order=T))
              col_color<-group_colorset(length(unique(col_color_index)))[col_color_index]
              rownames(expr_heatmap_data)<-expr_heatmap_data$Short_name
              expr_heatmap_data<-as.matrix(expr_heatmap_data[,ht_markers])
              #expr_heatmap_data<-apply(expr_heatmap_data,2,remove_extremum)
              expr_heatmap_data<-t(expr_heatmap_data)
              
              write.csv(expr_heatmap_data,paste0("./",output_dir,"/",cluster_name,"_csvs/",cluster_name,cluster_i,"_",key.title,file_tag,".csv"))
  
              
              x.inv<-try(heatmap.2(expr_heatmap_data,
                                   Rowv=Rowv,Colv=Colv,dendrogram=dendrogram,
                                   scale="none",
                                   key=T,keysize =keysize,key.title = key.title,
                                   trace="none",
                                   density.info="none",
                                   col=heatmap_color,
                                   srtCol=srtCol,
                                   srtRow=srtRow,
                                   margins =fig_margin,
                                   colCol = col_color,
                                   ColSideColors = col_color,
                                   main=paste0(key.title,"Heatmap",file_tag,"\n"),
                                   na.color=na.color),silent = T)
              
              if ('try-error' %in% class(x.inv)) {
                message("Warning- Could not output",key.title,"\n")
                paste(expr_heatmap_data)
              }
              
            }    
            
            pdf(file=paste0("./",output_dir,"/",cluster_name,'_',cluster_i,"_heatmap",file_tag,".pdf"),
                width=max(nrow(cluster_expr_trans[[cluster_i]])/5+4,8)*zoom_fig_width_h,
                height=max(length(ht_markers)/4+3,8)*zoom_fig_height_h
            )
            
            draw_heatmap(expr_trans,
                         heatmap_color=heatmap1_color,
                         key.title=paste0("Acsinh_Expr"))
            if(heatmap1_trans_method!="simpleAsinh") {
                draw_heatmap(expr_trans_ct,
                             heatmap_color=heatmap1_color,
                             key.title=paste0(heatmap1_trans_method,"_Expr"))}
            draw_heatmap(expr_trans_ratio,
                         heatmap_color=heatmap2_color,
                         key.title=paste0("Acsinh_Ratio_Expr"))
            
            
            dev.off()
            
        #======================各组的细胞数量统计=====================================================
            mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                             legend.key = element_rect(fill = "white", colour = "white"), #图标
                             legend.background = (element_rect(colour= "white", fill = "white")))
        
            cell_count<-cluster_cell_count[[cluster_i]] %>%
                        dplyr::filter_at(vars(one_of(major_cond)),all_vars(.%in% groups_to_show))
            cell_count$Short_name<-factor(cell_count$Short_name,levels=cell_count$Short_name,ordered=TRUE)
            #cell_count$n<-cell_count$n+1
        

            
            if(!exists("scale_y_limit_bar"))  scale_y_limit_bar=c(1,10*max(cell_count$n))
            if(!exists("scale_y_breakstep_bar"))   scale_y_breakstep_bar=50
            
            count_map<-ggplot(cell_count,aes_string(x="Short_name",y="n",fill=color_cond))+
                      geom_bar(stat = "identity")+
                      geom_hline(yintercept = 51, linetype=hline_line_type_bar,color = hline_line_color_bar, size=hline_line_size_bar)+
                      labs(y="number of cells Log(10)",title=paste0("Cell number of population", cluster_i," from each sample"))+
                      theme(axis.text= element_text(angle=45,hjust =1,vjust = 0.5,color = col_color))+
                      scale_y_log10(limits=scale_y_limit_bar,minor_breaks =scale_y_breakstep_bar)+
                      #scale_fill_brewer(palette = "Set1")
                      scale_fill_manual(values=brewer_color_sets(length(unique(col_color_index))))
            
            eval(parse(text = paste0("count_map<-count_map+",paralist_bar$themes)))
            
            pdf(file=paste0("./",output_dir,"/",cluster_name,"_",cluster_i,"_cell_number_barplot",file_tag,".pdf"),
                width=max(nrow(cluster_expr_trans[[cluster_i]])/5+4,8)*zoom_fig_width_bar,
                height=max(length(ht_markers)/4+3,8)*zoom_fig_height_bar
                )
            
            multiplot(count_map)
            dev.off()
            
            
        
        #======================绘制boxplot===================================================================
            
            #计算boxplot图片尺寸
            if(is.null(comparisons)) comparisons_n=2 else
              comparisons_n<-length(comparisons)+2
        
            boxplot_n<-length(ht_markers)
            major_cod_ID<-colnames(expr_trans)==major_cond
            cond_n<-length(unique(expr_trans[,major_cod_ID,drop=T]))
            singlewidth<-cond_n*0.5+1
            singleheight<-4+(comparisons_n-2)*0.2
        
            boxplot_ncol<-ceiling(sqrt(boxplot_n*4*singleheight/3/singlewidth))
            boxplot_nrow<-ceiling(boxplot_n/boxplot_ncol)
        
            boxplot_width<-boxplot_ncol*singlewidth
            boxplot_height<-boxplot_nrow*singleheight
        
            #各个marker的expression boxplot
            
            if(is.null(color_cond)) color_cond<-major_cond
            ncolor<-nrow(unique(expr_trans[,color_cond]))
            
            #定义boxplot函数
            

            if(!exists("border_size_b")) border_size_b=1
            if(!exists("border_color_b")) border_color_b="black"
            if(!exists("relative_y_limit_b")) relative_y_limit_b=c(0.1,0.1)
            
            marker_boxplot<-function(boxplot_marker,boxplot_data,ylab){
             
                          ydata<-boxplot_data[,boxplot_marker]
                          yrange<-max(ydata)-min(ydata)
                    
                          boxplot_p<-ggplot(boxplot_data,aes_string(x=major_cond,y=boxplot_marker,fill=color_cond))+
                            geom_boxplot(outlier.shape= NA,size=border_size_b,color=border_color_b,width=box_width_b)+
                            #geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
                            geom_jitter(shape=16, position=position_jitter(0.2))+
                            mytheme+
                            labs(title=boxplot_marker)+
                            scale_y_continuous(limits=c((min(ydata)-relative_y_limit_b[1]*yrange),(max(ydata)+relative_y_limit_b[2]*yrange*comparisons_n)))+
                            #scale_fill_brewer(palette = "Set1")+
                            scale_fill_manual(values=group_colorset(ncolor))+
                            labs(y=ylab)+
                            theme(legend.position = "none")+
                            #theme(axis.text= element_text(angle=45,hjust =1,vjust = 0.5,size=7))+
                            theme(axis.title = element_text(size=8))+
                            theme(title =element_text(size=8,face="bold"))
                          
                          eval(parse(text = paste0("boxplot_p<-boxplot_p+",paralist_boxplot$themes)))
                          
                          if(show_pvalues & is.null(comparisons)){
                            boxplot_p<-boxplot_p+
                              stat_compare_means(method = stat.method,paired=stat.paired,ref.group=heatmap_ctrl,aes(label = paste0("p = ", ..p.format..)))}
                          if(show_pvalues & !is.null(comparisons)){
                            boxplot_p<-boxplot_p+
                              stat_compare_means(method = comparisons.stat.method,paired=comparisons.stat.paired,comparisons=comparisons,aes(label = paste0("p = ", ..p.format..)))}
                          return(boxplot_p)
                        }
        
            all_boxplot_list<-lapply(ht_markers,marker_boxplot,boxplot_data=expr_trans,ylab="Expression Level")
            all_ratio_boxplot_list<-lapply(ht_markers,marker_boxplot,boxplot_data=expr_trans_ratio,ylab="Expression Arcsih Ratio to Ctrl")
            
            
            

                pdf(file=paste0("./",output_dir,"/",cluster_name,"_",cluster_i,"_Boxplots",file_tag,".pdf"),
                #pdf(file=paste0("./",output_dir,"/",cluster_name,cluster_i,"_",Manual_name,"(Symbol_", shape_char[cluster_i],") Boxplots",file_tag,".pdf"),
                width=boxplot_width*zoom_fig_width_b,
                height=boxplot_height*zoom_fig_height_b)
                multiplot(plotlist=all_boxplot_list,cols = boxplot_ncol)
                multiplot(plotlist=all_ratio_boxplot_list,cols = boxplot_ncol)
                
            dev.off()
        #==========================
    }
    
    if(length(heatmap_ctrl)==1){hide_ctrl=F}
    
    if(is.null(hide_ctrl)){
      
      if(Colv==T){
           message(paste0("Generating figures with and without heatmap_ctrl groups\n"))
           expr_report(hide_ctrl = T)
            expr_report(hide_ctrl = F)
      }else{
            expr_report(hide_ctrl = F) 
      }     
    }else{
      if(hide_ctrl==T)  {
            message(paste0("hide_ctrl is set to TRUE,heatmap_ctrl groups were hided in all output figures\n"))}
            expr_report(hide_ctrl)
    }

  }
  return_expr_list<-list()
  return_expr_list[["expr_ctrl_merge"]]<-expr_ctrl_merge
  return_expr_list[["expr_trans_merge"]]<-expr_trans_merge
  return(return_expr_list)
}



#'@title 生成marker表达的火山图
#'@description
#'依据stat_by_cluster生成的list数据，生成火山图，图中颜色代表不同的Marker，数字和字母等符号代表不同的cluster；
#'
#'@param cluster_stat       stat_by_cluster()函数生成的list数据
#'@param cond1              major_cond的一个，火山图（Volcano map）靠左一侧的条件
#'@param cond2              major_cond的一个，火山图（Volcano map）靠右一侧的条件
#'@param Usecolorset        火山图采用的配色方案 默认dif_seq_rainbow, 其他取值：rainbow，primary.color或者其他函数
#'@param show_adjust_p      一个逻辑变量， 取值TRUE时火山图显示校正后的p值，取值FALSE时显示原始p值
#'@param cluster_id         指定输出cluster的id，例如：如果要输出前五个cluster，cluster_id=c(1:5); 如果要输出全部cluster，保持其默认值NULL即可。
#'@param xlimit             一个包含两个数字的向量（vector），例如c(-5,+5),用于手动指定火山图横坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param ylimit             一个包含两个数字的向量（vector），例如c(0,3.5),用于手动指定火山图纵坐标范围，默认值为NULL，此时会自动设置xlimit
#'@param dif.level          一个数字，手动设置两组marker差异的的阈值大小，例如取默认值2时，两组均值相差2倍以上算显著（简单的说就是设置火山图上的两条竖线位置）
#'@param output_dir         输出数据文件夹名称
#'@param zoomfactor         一个正数，默认为1，控制图片的整体放缩比例 
#'@param label_repel        逻辑变量，label相互排斥分布，防止重叠  
#@param fig1_size,fig2_size,fig3_size  vector,规定图片的宽度和高度，例如c(10,5),宽度是10，高度是5 (单位：inch),与zoomfactor一起输出图片大小. 
#'                          fig1_size对应volcano，fig2_size对应volcano_by_markers, fig3_size 对应volcano_by_cluster

#'如使用全局统计默认参数，以下参数无需设置
#' @param stat.paired     逻辑变量，TRUE或者FALSE，表示是否是成对检验
#' @param conf.level      显著性水平，默认0.95,（即p<0.05为显著），调整该数值可以调整火山图上的横线位置
#' @param stat.method     一个字符串，指定统计分析的方法，取值有三种："auto"，“t-test"，“wilcox-test”,当取值为"auto"时，即为测试模式，会对数据进行正态分布和方差齐性测试，自动选择适合的统计方法，为正式分析过程中的统计方法选择提供依据。
#' @param set.equal.var   一个字符串，指定方差齐性，取值有三种："TRUE","FALSE"(默认),"auto",选择"auto"时即可以进行方差齐性测试；
#' @param p.adjust.method 对p值进行校正的方法："BH"(默认，推荐) "holm", "hochberg", "hommel", "bonferroni","BY","fdr", "none"
#' @return                none
#' @export


draw_expr_volcano<-function(cluster_stat,
                            cond1,
                            cond2,
                            Usecolorset=dif_seq_rainbow,
                            color_by_significant=FALSE,
                            show_adjust_p=FALSE,
                            cluster_id=NULL,
                            xlimit=NULL,
                            ylimit=NULL,
                            dif.level=2,
                            output_dir=paste0("Expr_volcano_plots_",cond1,"_vs_",cond2),
                            para_dir="0_Parameters",
                            #zoom_factor=1,
                            label_repel=F,
                            #  =c(7,8),
                            # fig2_size=c(20,20),
                            # fig3_size=c(20,20),

                            #如使用全局统计默认参数，以下参数无需设置
                            stat.paired=NULL,
                            stat.method=NULL,
                            set.equal.var=NULL,
                            conf.level=NULL,
                            p.adjust.method=NULL
                            ){


  
  if(0){
    
    Usecolorset=dif_seq_rainbow
    color_by_significant=FALSE
    show_adjust_p=FALSE
    cluster_id=NULL
    xlimit=NULL
    ylimit=NULL
    dif.level=2
    output_dir=paste0("Expr_volcano_plots_",cond1,"_vs_",cond2)
    para_dir="0_Parameters"
    zoom_factor=1
    fig1_size=c(15,8)
    fig2_size=c(20,20)
    fig3_size=c(20,20)
    
    #如使用全局统计默认参数，以下参数无需设置
    stat.paired=NULL
    stat.method=NULL
    set.equal.var=NULL
    conf.level=NULL
    p.adjust.method=NULL
    
    }
 

  if (!dir.exists(paste0("./",output_dir))) {
    dir.create(paste0("./",output_dir))
  }

  copy_para_csv(output_dir,c("expr_volcano_barplot.csv","expr_volcano_volcanoplot.csv","expr_volcano.pptx"))
 
  wdir=getwd()  
  alias_as_label               <- cluster_stat[["metadata"]]$alias_as_label
  clustern               <- cluster_stat[["metadata"]]$clustern
  summerise_method       <- cluster_stat[["metadata"]]$summerise_method
  major_cond             <- cluster_stat[["metadata"]]$major_cond
  all_markers            <- cluster_stat[["all_markers"]]
  merged_cluster_id      <-cluster_stat[["merged_cluster_id"]]



  cat(paste0("Output to folder: ",wdir,"/",output_dir,"\n"))

  ht_markers<-as.character(subset(all_markers,expr_para==1)$markers)

  major_cond<-as.character(major_cond)
  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)
  

  cluster_expr_raw  <-cluster_stat[["cluster_expr_raw"]]
  if(is.null(stat.paired)) stat.paired        <- cluster_stat[["metadata"]]$stat.paired
  if(is.null(stat.method)) stat.method        <- as.character(cluster_stat[["metadata"]]$stat.method)
  if(is.null(conf.level)) conf.level          <- cluster_stat[["metadata"]]$conf.level
  if(is.null(set.equal.var)) set.equal.var    <- as.character(cluster_stat[["metadata"]]$set.equal.var)
  if(is.null(p.adjust.method)) p.adjust.method    <- as.character(cluster_stat[["metadata"]]$p.adjust.method)

  if(is.null(cluster_id)){
    cluster_id=merged_cluster_id
  }

  if(summerise_method=="mean"){
    use_method<-mean}else
      if(summerise_method=="median")
      { use_method<-median
      }


  if(show_adjust_p==TRUE){
    p_para<-"log_p_adjust"
    adjust_label="Adjust "}else{
    p_para<-"log_p_values"
    adjust_label=""
    }


  if(stat.method=="auto") stat.method="t-test or wilcox-test\n(Statistical method is selected automaticly)"

  if(stat.paired){
    if_paired<-"paired"}else{
      if_paired<-"unpaired"}


  cond_name<-major_cond
  cond_name<-as.character(cond_name)
  #画火山图
  volcano_data_all<-list()

  for(cluster_i in cluster_id){

    expr_raw<-cluster_expr_raw[[cluster_i]]

  #使用Acsinh transformed values 计算Pvalue
    #unpaired_control_cond<-dplyr::filter_at(expr_trans_cond,vars(one_of(cond_name)),all_vars(.==heatmap_ctrl_i))
    
  cond1_data<-
    expr_raw %>%
    #dplyr::filter_(paste0(cond_name,"==cond1"))%>%
    dplyr::filter_at(vars(one_of(cond_name)),all_vars(.==cond1))%>%
    dplyr::select(ht_markers)
  cond2_data<-
    expr_raw %>%
    #dplyr::filter_(paste0(cond_name,"==cond2"))%>%
    dplyr::filter_at(vars(one_of(cond_name)),all_vars(.==cond2))%>%
    dplyr::select(ht_markers)


  volcano_cluster_data<-combined_stat(cond1_data,
                                      cond2_data,
                                      alternative="two.sided",
                                      mu=0,
                                      transform.method="asinh",
                                      stat.paired=stat.paired,
                                      conf.level=conf.level,
                                      stat.method=stat.method,
                                      set.equal.var=set.equal.var,
                                      p.adjust.method=p.adjust.method,
                                      silent=TRUE)


  #计算Ratio

  cond1_data_m<-
    cond1_data %>%
    summarise_if(is.numeric,use_method,na.rm=F)

  #summarise_each(funs(use_method))

  cond2_data_m<-
    cond2_data %>%
    summarise_if(is.numeric,use_method,na.rm=F)
  #summarise_each(funs(use_method))

  No_zero<-0.2  #add to prevent divided by zero

  log_ratio=log2((cond2_data_m+No_zero)/(cond1_data_m+No_zero))



  volcano_data<-data.frame(marker=colnames(cond1_data),
                           cluster_i=cluster_i,
                           cluster_shape=shape_char[cluster_i],
                           volcano_cluster_data,
                           log_ratio=t(log_ratio)
                          )

  if(cluster_i==1){volcano_data_all=volcano_data}else{
    volcano_data_all<-rbind(volcano_data_all,volcano_data)
  }

}
  volcano_data_all$log_p_values=-log10(volcano_data_all$p_values)
  volcano_data_all$log_p_adjust=-log10(volcano_data_all$p_adjust)

  # if(!is.null(Manual_Pop_Index)){
  #       volcano_data_all$cluster<-paste0("cluster",volcano_data_all$cluster_i)
  #       volcano_data_all<-dplyr::left_join(volcano_data_all,Manual_Pop_Index,by="cluster")
  #       volcano_data_all$cluster<-NULL
  #       }
  cat(paste0("Output Statistical Report.csv...\n"))
  write.csv(volcano_data_all,paste0("./",output_dir,"/","Statistical Report.csv"),row.names = TRUE)

  shape_index<-c(49:57,65:90,97:122)  #数字，大写字母，小写字母  shape
  shape_char<-intToUtf8(shape_index,multiple=T)

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))


  #准备FDR-pValue画图
  FDR_ranks<-factor(c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),levels=c("other(adjust_p>=10%)","adjust_p<10%","adjust_p<5%"),ordered=TRUE)
  FDR_pValue<-dplyr::arrange(volcano_data_all,p_values)
  FDR_pValue$row_id<-factor(c(1:nrow(FDR_pValue)))
  FDR_pValue$FDR_rank<-FDR_ranks[1]
  FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.1]=FDR_ranks[2]
  FDR_pValue$FDR_rank[FDR_pValue$p_adjust<0.05]=FDR_ranks[3]
  cat(paste0("Output Statistical Report(arranged).csv...\n"))
  write.csv(FDR_pValue,paste0("./",output_dir,"/","Statistical Report(arranged).csv"),row.names = TRUE)

  FDR_table_data<-FDR_pValue[(FDR_pValue$p_values<0.05)&(!is.na(FDR_pValue$p_values)),c(1:6,8,10)]
  round_2<-function (var){return(round(var,3))}
  if(nrow(FDR_table_data)>0)
  FDR_table_data[,c(4,5,7)]<-apply(FDR_table_data[,c(4,5,7)],2,round_2)

  FDR_table_plot<-ggtable(FDR_table_data,show_row_name = F)


  FDR_pValue<-dplyr::arrange(FDR_pValue,desc(p_values))
  FDR_pValue$row_id<-factor(c(1:nrow(FDR_pValue)))
  FDR_pValue$p_values=-log10(FDR_pValue$p_values)
# FDR_pValue$FDR_rank<-factor(FDR_pValue$FDR_rank)

  #FDR_pValue<-na.omit(FDR_pValue)
  FDR_pValue <-FDR_pValue[complete.cases(FDR_pValue[,c("p_values","p_adjust","log_ratio","log_p_values")]),]
  
  paralist2<-get_para(output_dir,"expr_volcano_barplot.csv")

  #paralist2<-get_para(para_file = paste0("./",output_dir,"/para_csvs/expr_volcano2.csv"))
  eval(parse(text =paralist2$items))
  if(!exists("scale_y_breakstep")) scale_y_breakstep2<-max(0.5,ceiling(max(FDR_pValue$p_values/5)))
  if(!exists("ylimit2")) ylimit2<-c(0,ceiling(max(FDR_pValue$p_values)))

  FDR_pValue_plot<-ggplot(FDR_pValue)+
    geom_bar(aes_string(x="row_id",y="p_values",fill="FDR_rank"),stat="identity",width =bar_width2)+
    labs(y="-Log10(p_values)",title="Summary of Adjust p Value")+
    labs(x=paste0(""))+
    coord_flip()+
    labs(subtitle=paste0("Adjust p values were separately calculated based on the p values of each cluster "))+
    mytheme+
    #expand_limits(y=c(0,3))+
    theme(panel.background=element_rect(colour='white',size = 1 ))+
    theme(axis.line.x=element_line(colour='black',linetype = NULL,size=0.3))+
    theme(axis.ticks.y=element_blank())+
    scale_fill_manual(values = bar_color2)+
    scale_y_continuous(breaks=seq(ylimit2[1],ylimit2[2],by=scale_y_breakstep2),minor_breaks = NULL)+
    theme(legend.position = c(0.7,0.4),legend.title=element_blank())+
    theme(axis.text.y=element_blank())+
    theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))
  eval(parse(text = paste0("FDR_pValue_plot<-FDR_pValue_plot+",paralist2$themes)))
  

  #为了方便绘图，所有小于0.0001的p值等于0.0001
  if(is.null(xlimit)){
    xuplimit=max(abs(volcano_data_all$log_ratio)+2,na.rm = TRUE)
    xlimit=c(-1*xuplimit,xuplimit)}
  
  
  Y_ID<-colnames(volcano_data_all)==p_para
  
  if(is.null(ylimit)){
    yuplimit=min(max(volcano_data_all[Y_ID]+1.5,na.rm = TRUE),4)
    ylimit=c(0,yuplimit)}  
  
  yuplimit<-ylimit[2]
  volcano_data_all$marker<-as.character( volcano_data_all$marker)
  volcano_data_all$p_values[volcano_data_all$p_values<10^(-yuplimit)]=10^(-yuplimit)
  volcano_data_all$p_adjust[volcano_data_all$p_adjust<10^(-yuplimit)]=10^(-yuplimit)
  volcano_data_all$log_p_values=-log10(volcano_data_all$p_values)
  volcano_data_all$log_p_adjust=-log10(volcano_data_all$p_adjust)
  volcano_data_all$significant<-TRUE
  significant_id<-unique(c(which(volcano_data_all$p_values>(1-conf.level)),
                           which(abs(volcano_data_all$log_ratio)<1),
                           which(is.na(volcano_data_all$p_values))))
  volcano_data_all[significant_id,"significant"]<-FALSE


  shape_list<-data.frame(cluster_i=merged_cluster_id,
                         shape_index=shape_index[merged_cluster_id],
                         cluster=paste0("cluster",merged_cluster_id)) 
  if(alias_as_label){
         volcano_data_all$cluster_i<-factor(volcano_data_all$cluster_i,levels = levels(merged_cluster_id),ordered = T)}
  volcano_data_all<-full_join(volcano_data_all,shape_list,by="cluster_i")
  # if(!is.null(Manual_Pop_Index)){
  #   
  #   shape_list<-left_join(shape_list,Manual_Pop_Index,by="cluster")
  #   
  # }

  
  #设置色彩：
  marker_n<-length(unique(volcano_data_all$marker))
  
  color_ref<-data.frame(marker=unique(volcano_data_all$marker),
                        raw_colors=Usecolorset(marker_n),
                        stringsAsFactors = F
                        )
  volcano_data_all<-left_join(volcano_data_all,color_ref,by="marker")
  volcano_data_all$sig_colors<-volcano_data_all$raw_colors

    for (this_row in row.names(volcano_data_all)) {
    if(!volcano_data_all[this_row,"significant"]){
                          volcano_data_all[this_row,"sig_colors"]<-ColToGray(volcano_data_all[this_row,"raw_colors"])
                          # volcano_data_all[this_row,"sig_colors"]<-colorRampPalette(colors=c(volcano_data_all[this_row,"raw_colors"],"white"))(4)[3]#过度色
                         }
    }
                         
  raw_color_list<-unique(volcano_data_all$raw_colors)
  #volcano_data_all$raw_colors<-factor(volcano_data_all$raw_colors,levels=raw_color_list,ordered=T)
  marker_list<-unique(volcano_data_all$marker)
  volcano_data_all$marker<-factor(volcano_data_all$marker,levels=marker_list,ordered=T)
  
  
  sig_color_list<-unique(volcano_data_all$sig_colors)
  volcano_data_all$sig_colors<-factor(volcano_data_all$sig_colors,levels=sig_color_list,ordered=T)

  mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                   legend.key = element_rect(fill = "white", colour = "white"), #图标
                   legend.background = (element_rect(colour= "white", fill = "white")))


   cat(paste0("Output Volcanoplot...\n"))

  x_threshold<-round(log2(dif.level),digits=1)

  paralist1<-get_para(output_dir,"expr_volcano_volcanoplot.csv")
  eval(parse(text =paralist1$items))
  if(!exists("xlimit1")) xlimit1<-xlimit
  if(!exists("ylimit1")) ylimit1<-ylimit
  
  volcanoplot<-ggplot(volcano_data_all)+
        #scale_shape_identity() +
        #geom_point(aes_string(x="log_ratio",y=p_para,colour="marker",shape="shape_index"),size=4)+
        geom_point(aes_string(x="log_ratio",y=p_para,colour="marker"),size=dotsize1)+
        labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
        labs(x=paste0("Log2_Ratio","\n", cond1, " <- versus -> ", cond2))+
        labs(subtitle=paste0(adjust_label,"p values of arcsinh transformed expression were calculated by ",if_paired," ",stat.method))+
        mytheme+
        expand_limits(x=xlimit,y=ylimit)+
        geom_vline(xintercept = c(x_threshold,-1*x_threshold), linetype=vline_line_type1,color = vline_line_color1, size=vline_line_size1)+
        geom_hline(yintercept=-log10(1-conf.level), linetype=hline_line_type1, color =hline_line_color1,size=hline_line_size1)+
        scale_x_continuous(breaks=c(seq(floor(xlimit1[1]),ceiling(xlimit1[2]),by=1)),minor_breaks = NULL)+
        scale_y_continuous(breaks=c(seq(floor(ylimit1[1]),ceiling(ylimit1[2]),by=1),1.3),minor_breaks = NULL)+
        theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
        theme(panel.background=element_rect(colour='black',size = 1 ))+
        theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
        scale_color_manual(values = raw_color_list)
  if(!label_repel){
        volcanoplot<-volcanoplot+geom_text(aes_string(x="log_ratio",y=p_para,label="cluster_i"),colour=label_color1,size=lable_size1 )
  }else{
        volcanoplot<-volcanoplot+geom_text_repel(aes_string(x="log_ratio",y=p_para,label="cluster_i"),colour=label_color1,size=lable_size1 )
  }
    eval(parse(text = paste0("volcanoplot<-volcanoplot+",paralist1$themes)))
    
  volcanoplot_legend<-as_ggplot(get_legend(volcanoplot))
  volcanoplot<-volcanoplot+theme(legend.position = "none")
  
if(color_by_significant){
    
  volcanoplot<-ggplot(volcano_data_all)+
        #scale_shape_identity() +
        #geom_point(aes_string(x="log_ratio",y=p_para,colour="marker",shape="shape_index"),size=4)+
        geom_point(aes_string(x="log_ratio",y=p_para,colour="sig_colors"),size=dotsize1)+
        labs(y="-Log10(p_values)",title="Differential Expression of Markers")+
        labs(x=paste0("Log2_Ratio","\n", cond1, " <- versus -> ", cond2))+
        labs(subtitle=paste0(adjust_label,"p values of arcsinh transformed expression were calculated by ",if_paired," ",stat.method))+
        mytheme+
        expand_limits(x=xlimit,y=ylimit)+
        geom_vline(xintercept = c(x_threshold,-1*x_threshold), linetype=vline_line_type1,color = vline_line_color1, size=vline_line_size1)+
        geom_hline(yintercept=-log10(1-conf.level), linetype=hline_line_type1, color =hline_line_color1,size=hline_line_size1)+
        scale_x_continuous(breaks=c(seq(floor(xlimit1[1]),ceiling(xlimit1[2]),by=1)),minor_breaks = NULL)+
        scale_y_continuous(breaks=c(seq(floor(ylimit1[1]),ceiling(ylimit1[2]),by=1),1.3),minor_breaks = NULL)+
        theme(plot.margin=unit(c(1,1,1,1), unit = "cm"))+
        theme(panel.background=element_rect(colour='black',size = 1 ))+
        theme(axis.line=element_line(colour='black',linetype = 1,size=0.3))+
        scale_color_manual(values = sig_color_list)+
#        geom_text(aes_string(x="log_ratio",y=p_para,label="cluster_i"),colour=label_color1,size=lable_size1 )+
        theme(legend.position="none")
   if(!label_repel){
        volcanoplot<-volcanoplot+geom_text(aes_string(x="log_ratio",y=p_para,label="cluster_i"),colour=label_color1,size=lable_size1 )
  }else{
        volcanoplot<-volcanoplot+geom_text_repel(aes_string(x="log_ratio",y=p_para,label="cluster_i"),colour=label_color1,size=lable_size1 )
  }
  
   eval(parse(text = paste0("volcanoplot<-volcanoplot+",paralist1$themes)))
  
  
}
       
  volcano_list<-list(volcano_data_all,volcanoplot)
  # save(volcano_list,
  #     file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,"_1.RData"))
  
# if(!is.null(Manual_Pop_Index)){
#       
#         shape_table<- data.frame(label=merged_cluster_id,
#                                  Population=shape_list$Pop_name)
#         volcano_note<-ggtable(shape_table,show_row_name=F,show_col_name=T)+
#                       theme(plot.margin=unit(c(0,0,0,0), unit = "cm"))
#         
#         layout<-cbind(matrix(c(0,rep(1,4),rep(1,4),0),ncol = 1, byrow = T),
#                       matrix(c(rep(2,40)),ncol = 4, byrow = T),
#                       matrix(c(0,rep(3,4),rep(3,4),0),ncol = 1, byrow = T),
#                       matrix(c(rep(4,40)),ncol = 4, byrow = T))
#         
#         pdf(file=paste0("./",output_dir,"/","Volcano_",cond1," vs ",cond2,".pdf"), width=zoom_factor*fig1_size[1], height=zoom_factor*fig1_size[2])
#         
#         multiplot(volcano_note,volcanoplot,volcanoplot_legend,FDR_pValue_plot,layout=layout)
#         
#         dev.off()
# }else{
  # layout<-cbind(matrix(c(rep(1,40)),ncol = 4, byrow = T),
  #               matrix(c(0,rep(2,4),rep(2,4),0),ncol = 1, byrow = T),
  #               matrix(c(rep(3,40)),ncol = 4, byrow = T))
  # 
  
  pdf(file=paste0("./",output_dir,"/","1_Volcano_",cond1," vs ",cond2,".pdf"), width=zoom_fig_width11*8, height=zoom_fig_height11*10)
  
  volcanoplot_merge<-ggarrange(volcanoplot,volcanoplot_legend,#gblank,gblank,
                               ncol = 2,
                               widths = c(1,legend_wide1),
                               heights = 1)
  #multiplot(volcanoplot,volcanoplot_legend,FDR_pValue_plot,layout=layout)
  
  dev.off()
  
  
  pdf(file=paste0("./",output_dir,"/","4_barplot_",cond1," vs ",cond2,".pdf"), width=zoom_fig_width2*6, height=zoom_fig_height2*10)

      multiplot(FDR_pValue_plot)
  
  dev.off()

#     }

  
  cat(paste0("Output Volcanoplot facet by markers...\n"))

  pdf(file=paste0("./",output_dir,"/","2_Volcano_",cond1," vs ",cond2,"by markers.pdf"),
      width=zoom_fig_width21*20, height=zoom_fig_height21*20)

      volcanoplot2<-volcanoplot+facet_wrap(~marker)
  
  multiplot(volcanoplot2)
  dev.off()


  cat(paste0("Output Volcanoplot facet by clusters...\n"))

    # if(!is.null(Manual_Pop_Index)){
    #   colour_channel<-"Pop_name"
    # }else{
      colour_channel<-"cluster_i"
 #   }

 
    # if(!is.null(Manual_Pop_Index)){
    #   volcanoplot3<-volcanoplot+facet_wrap(~Pop_name)
    #   
    # }else{
      volcanoplot3<-volcanoplot+facet_wrap(~cluster_i)
    #}
    
  pdf(file=paste0("./",output_dir,"/","3_Volcano_",cond1," vs ",cond2,"by population.pdf"),
      width=zoom_fig_width31*20, height=zoom_fig_height31*20)
  
     multiplot(volcanoplot3)
  dev.off()

}



#' @title Arcsinh 转化函数

#' @description：
#' 这是一个数据转化函数，由Cytofkit::cytofAsinh() 函数修改而来，根据需要进行了简化，主要去掉了原函数关于扣除背景和随机化的部分。
#'
#' @param value the data need to be transformed,a vector of numeric values.
#' @param cofactor Cofactor for asinh transformation, default 5 for mass cytometry data.
#' @return a vector of numeric values containing transformed data
#' @export

simpleAsinh <- function(value, cofactor = 5) {
  value <- value / cofactor
  value <- asinh(value)
  return(value)
}



#' @title 去除极值
#' @description：
#' 去除极值，将大于0.95分位数的数值赋值为0.95分位数，小于0.05分为数的数值赋值为0.05分位数
#' @param value 待处理的数据
#' @return a vector of numeric values containing treated data
#' @export
#'
remove_extremum<-function(value,uppct=0.95,lowpct=0.05){
  quantile_95<-quantile(value,uppct,na.rm=TRUE)
  quantile_05<-quantile(value,lowpct,na.rm=TRUE)
  value[value>quantile_95]<-quantile_95
  value[value<quantile_05]<-quantile_05
  return(value)
}

