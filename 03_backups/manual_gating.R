
#' @title generate_gating_ref
#'
#' @description
#' 生成每个文件的manual gating ref文件.   ref文件为csv文件，其中行为细胞，列为每个gate，1代表细胞在这个gate中，0代表不在。

#'
#' @param gs               包含gating信息和单细胞表达数据的gateset object
#' @param fig_output_dir   函数会导出每个FCS文件的Gating strategy和Gating hearacy图片,figoutputdir是图片所在的子文件夹名称 
#' @param csv_output_dir   gating_ref所在的目录，默认和原始FCS文件放在一起
#' @param rawfcsdir        原始FCS文件所在目录 
#' 
#' 
#' @return                 无返回值
#' @export

generate_gating_ref<-function(  gs,
                                fig_output_dir="gating hierarcy",
                                csv_output_dir="downsample",
                                rawfcsdir=raw_fcs_dir){
  
  if(0){
    
    fig_output_dir="gating hierarcy"
    csv_output_dir=raw_fcs_dir
    
  }
  
  if (!dir.exists(paste0("./",csv_output_dir))) {
    dir.create(paste0("./",csv_output_dir))
  }
  
  
  if (!dir.exists(paste0("./",fig_output_dir))) {
    dir.create(paste0("./",fig_output_dir))
  }
  
  file_name <- list.files(rawfcsdir,pattern='.fcs$', full=TRUE)
  File_Id<-sub("^.*/","",file_name)
  File_Id<-sub(".fcs$","",File_Id)
  
  fcs_num<-length(file_name)
  #gate_names<-getNodes(gs,path="auto")
  gate_names<-gs_get_pop_paths(gs,path = "auto")
  ngates<-length(gate_names)
  figdim<-ceiling(sqrt(ngates))
  filter_table_list<-list()
  #cat(paste0("Outputing the : ",rawfcsdir,"\n"))
  
  
  for (fcs_id in c(1:fcs_num)) {
    
    message("Start to treatwith ",File_Id[fcs_id],":\n")
    if(file.exists(paste0(csv_output_dir,"/",File_Id[fcs_id],".csv"))){
         message(paste0("Warning: ",File_Id[fcs_id],".csv"," already exist,stop to avoid overide.\n"))
         next
    }
    
    this_ff<-gs@data[[fcs_id]]
    
    #this_ff <- read.FCS(file_name[fcs_id],transformation = FALSE,truncate_max_range=F)
    # asinhTrans <- arcsinhTransform(a=0, b=0.2, c=0)
    # translist <- transformList(names(gs[fcs_id]@transformation[[1]]), asinhTrans)
    # this_ff <- transform(this_ff, translist)
    filter_table<-as.data.frame(matrix(data = NA,nrow=nrow(this_ff@exprs),ncol=0))
    for (this_gate_name in gate_names[-1]) {
      #this_gate_name=gate_names[3]
      this_gate<-getGate(gs[[fcs_id]],this_gate_name)  
      fiter_result = flowCore::filter(this_ff,this_gate)
      filter_table[,this_gate_name]<-as.numeric(fiter_result@subSet)
    }
    
    #cat(apply(filter_table, 2, sum))
    # filter_table_list[[fcs_id]]<-filter_table
    write.csv(filter_table,
              file = paste0(csv_output_dir,"/",File_Id[fcs_id],".csv"),
              row.names = FALSE)
    cat(paste0("Outputed: ",csv_output_dir,"/",File_Id[fcs_id],".csv","\n"))
    tiff(file=paste0("./",fig_output_dir,"/",this_ff@description$GUID," gating hierarchy.tiff"),width =figdim*3*72 ,height = figdim*3*72 )
    a<-suppressMessages(autoplot(gs[[fcs_id]]))
    multiplot(a)
    #dev.off()
    graphics.off()
    cat(paste0("Outputed: ",fig_output_dir,"/",this_ff@description$GUID," gating hierarchy.tiff","\n\n"))
    
  }    
  cat(paste0("Outputing Gating Hierachy Plot...\n"))
  tiff(file=paste0("./",fig_output_dir,"/",gs@data[[1]]@description$GUID," gating hierarchy.tiff"),width =figdim*3*72 ,height = figdim*3*72 )
  plot(gs)
  #multiplot(hierarchy_plot)
  graphics.off()
  cat(paste0("Finished.\n"))
  
}


#' @title generate_gating_strategy
#'
#' @description
#' 生成一个csv文件，反映每个population与各个gate的关系。横坐标是population名称，纵坐标是Gate的名称。

#' @param gs               包含gating信息和单细胞表达数据的gateset object
#' @param csv_output_dir   gating_ref所在的目录，默认和原始FCS文件放在一起
#' 
#' 
#' @return                 无返回值
#' @export

#

generate_gating_strategy<-function(gs,
                                   csv_output_dir=metadata_dir){
  
  pop_table<-gs_pop_get_count_fast(gs[[1]])[,]
  pop_paths<-pop_table$Population
  pop_names<-sub("^.*/","",pop_paths)
  gating_strategy_table<-as.data.frame(matrix(data = NA,ncol = length(pop_paths)+1,nrow=0))
  colnames(gating_strategy_table)<-c("Pop_Extract_Seq",pop_names)
  
  
  for (pop_name_id in c(1:length(pop_names))) 
  {
    used_gates<-strsplit(pop_paths[pop_name_id],split = "/")
    used_gates<-used_gates[[1]]
    used_gates<-used_gates[-1]
    gating_strategy_table[pop_names[pop_name_id],-1]=0
    gating_strategy_table[pop_names[pop_name_id],used_gates]=1 
  }
  gating_strategy_table[,1]="" 
  
  if(file.exists(paste0(csv_output_dir,"/gating_strategy.csv"))){
    message("Warning: gating_strategy.csv already exist,stop to avoid overide.\n") 
  } else{
    write.csv(gating_strategy_table,
              file = paste0(csv_output_dir,"/gating_strategy.csv"))
    cat(paste0("Outputed: ",csv_output_dir,"/gating_strategy.csv","\n\n"))
    
  }
}







#' @title gating_ref_merge
#' @description
#' 根据提供的rowname reference, 抽取gating ref文件中的细胞信息，并合并

#'
#' @param gating_ref_csv    vector，包含gating_ref csv 文件的绝对目录
#' @param rowname_ref         vector，包含前面程序提取的所有细胞的行名，使用combined_data_raw的行名就行了

#' 
#' 
#' @return                 无返回值
#' @export
#' 
#' 
gating_ref_merge<-function(gating_ref_csv = gating_ref_csv#,
                           #  rowname=row.names(combined_data_raw)
                           ){
  
  combined_gating_ref<-read.csv(gating_ref_csv[1],header = T)
  combined_gating_ref<-combined_gating_ref[0,]
  for (this_gate_ref in gating_ref_csv) {
     cat(paste("Start to treat with ",sub("^.*/","",this_gate_ref),"...   "))
     #this_gate_ref=gating_ref_csv[1]
     gate_table<-read.csv(this_gate_ref,header = T)
     File_ID<-sub(".csv$","",basename(this_gate_ref))
     row.names(gate_table)<-paste0(File_ID,"_",row.names(gate_table))
     combined_gating_ref<-rbind(combined_gating_ref,
                                  gate_table)
     cat(paste("Finished.\n"))
     
     
  }
  
  #combined_gating_ref<-combined_gating_ref[rowname,]
  
  # if(sum(is.na(rownames(combined_gating_ref)))){
  #   combined_gating_ref<-NULL
  #   stop("Row names of csv and fcs is not compatible.\n")
  # }
  # 
  return(combined_gating_ref)
}






#' @title extract_pop_ref
#'
#' @description
#'  从combined_gating_ref中抽取单细胞的populatin信息，包含单细胞与每个population的从属关系和每个单细胞的manual gating population注释
#'  
#'  
#' @param combined_gating_ref   dataframe, gating_ref_merge产生的合并后的gating_ref数据
#' @param gating_strategy_table   dataframe，包含每个population与各个gate的关系以及所要提取的population顺序信息
#' @param output_dir              vector，导出Manual_Pop_Index.csv 的目录

 
#' @return                 无返回值
#' @export


extract_pop_ref<-function(combined_gating_ref,
                          gating_strategy_table,
                          output_dir="output"
                          ){

        combined_pop_ref<-combined_gating_ref[,0]
        
        gate2pop_fiter<-function(gate_vector){
                         if(all(gate_vector==T)) {return(1)}else{
                           return(0)
                         }
        }
        
        if(all(is.na(gating_strategy_table$Pop_Extract_Seq))){
          stop("Pop_Extract_Seq is not setted, please check the gating_strategy.csv.\n")
        }
        Pop_Extract<-arrange(gating_strategy_table[,c("X","Pop_Extract_Seq")],Pop_Extract_Seq)
        Pop_Extract<-Pop_Extract[!is.na(Pop_Extract$Pop_Extract_Seq),]

        layer_id<-which(colnames(gating_strategy_table)=="Pop_Extract_Seq")
        if(is.na(layer_id[1])) stop("Pop_Extract_Seq column is not detected in gating_strategy.csv\n")
        gating_strategy_table<-gating_strategy_table[,-layer_id]
        
        
        message(paste0("Start to treatwith selected populations:\n"))
        
                
        for(this_marker in gating_strategy_table[,1]){
              
              this_row<-which(gating_strategy_table[,1,drop=T]==this_marker)
              gating_cols<-which(gating_strategy_table[this_row,,drop=F]==1)
              adj_gating_cols<-gating_cols-1
              
              pop_ref<-apply(combined_gating_ref[,adj_gating_cols,drop=F],1,gate2pop_fiter)
              
              combined_pop_ref<-cbind(combined_pop_ref,pop_ref)}
        
        colnames(combined_pop_ref)<-gating_strategy_table[,1,drop=T]
        
        combined_pop_ref[,"manual_gate"]="Others"
        #combined_pop_ref[,"pop_index"]  = max(Pop_Extract$Pop_Extract_Seq)+1
        combined_pop_ref[,"pop_index"]=0
        for (this_pop in Pop_Extract$X) {
          for (thisrow in c(1:nrow(combined_pop_ref))) {
            if(combined_pop_ref[thisrow,this_pop]==1){
                combined_pop_ref[thisrow,"manual_gate"]<-this_pop
                combined_pop_ref[thisrow,"pop_index"]<-Pop_Extract[which(Pop_Extract$X==this_pop),"Pop_Extract_Seq"]}
              }
          cat(paste0("Population ",this_pop, " annotation finished.\n"))
        }
        message(paste0("Outputed:Extract_Pop_Index.csv.\n\n"))
        #Pop_Extract<-rbind(Pop_Extract,c("Others",max(Pop_Extract$Pop_Extract_Seq)+1))
        Pop_Extract<-rbind(Pop_Extract,c("Others",0))
        
        colnames(Pop_Extract)<-c("Pop_name","manual_gate")
        Pop_name<-Pop_Extract$Pop_name
        # Pop_Extract$Pop_name
        # 
        # Pop_name<-gsub("-","_",paraname)
        # Pop_name<-gsub("+","_",paraname)
        Pop_name<-gsub(" ","_",Pop_name)
        Pop_Extract$Pop_name<-Pop_name
        if (!dir.exists(paste0("./",output_dir))) {
          dir.create(paste0("./",output_dir))
        }
        
        
        message('Output cluster label list to following folder:\n ')
        cat(paste0("./",output_dir,"/manual_gate_alias.csv\n"))
        write.csv(Pop_Extract[,c(2,1)],
                  file = paste0("./",output_dir,"/manual_gate_alias.csv"),
                  row.names = F)
        message("All population annotation has been finished.\n\n")
        

        cat("Raw Manual Gated population statics:\n")
        raw_pop_stat<-apply(combined_pop_ref[,c(1:(ncol(combined_pop_ref)-2))],2, sum)
        cat(paste(names(raw_pop_stat), raw_pop_stat, sep = ":", collapse = "   "))
        
        cat("\n\nExtracted Manual Gated population statics:\n")
        
        Extracted_Pop_stat<-combined_pop_ref[,'manual_gate',drop=F]%>%
           group_by(manual_gate)%>%
             summarise(num=n())
        
        for (i in c(1:nrow(Extracted_Pop_stat))) {
          cat(paste0(Extracted_Pop_stat[i,1],":",Extracted_Pop_stat[i,2],"  "))
          
        }
        

        return(combined_pop_ref$pop_index)
}



#设置别名
generate_cluster_alias<-function(clustering_result=metacluster_result,
                                 cluster_name="metacluster"
                                 ){
  
                          if (!dir.exists(paste0("./output"))) {
                            dir.create(paste0("./output"))
                          }
                          cluster_alias_df<-data.frame(sort(unique(metacluster_result)),alias="")
                          colnames(cluster_alias_df)<-c(cluster_name,paste0(cluster_name,"_alias"))     
                          write.csv(cluster_alias_df,paste0("./output/",cluster_name,"_alias.csv"),row.names = F)
                          cat(paste0("./output/",cluster_name,"_alias.csv"," outputed.\n Please set the alias manually before use it."))
                          }
