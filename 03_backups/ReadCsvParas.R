

#'@title 获取para_csv 文件中的参数设置
#'@description
#'将para_csv中的参数设置读取出来，转化为可以运行的代码
#'
#' @param para_flie     字符串，csv文件的名称

#' @return              list变量，包含themes和items参数（themes就是ggplot2 themes的参数，其他参数在item里面）
#' @export
# para_path<-output_dir
# para_file<-"draw_tsne_figs_dotplots.csv"
# str(para_df)
# getsuffix("")

get_para<-function(para_path="",
                   para_file="",
                   theme_silent=T,
                   item_silent=T){
  
                   para_df<-read.csv(paste0("./",para_path,"/para_csvs/",para_file),
                                     header = T,stringsAsFactors = F)

                   
                   pvid<-which(colnames(para_df) %in%c("Parameters","Value","suffix"))
                   
                   theme_settings<-para_df[which(para_df$Value!=""& para_df$Type=="theme"),pvid]
                   
                   
                   getsuffix<-function(suffix){
                   if(is.na(suffix))suffix=""
                   if(is.null(suffix))suffix=""
                   return(suffix)
                   }
                   
                   
                   paras<-c()
                   for (theme_row in row.names(theme_settings)) {
                     
                     parai<-paste0(theme_settings[theme_row,"Parameters"],
                            "=",
                            theme_settings[theme_row,"Value"])
                   
                     paras<-c(paras,parai)
                     }
                   
                   
                   themes<-paste0("theme(",paras,")",collapse = "+")
                   
                   if(theme_silent==F) {cat(themes)}
                   para_df[para_df$Type=="function","suffix"]=""
                   item_settings<-para_df[which(para_df$Value!=""& para_df$Type!="theme"),pvid]
                   
                   paras2<-c()
                   for (item_row in row.names(item_settings)) {
                     
                     
                     parai2<-paste0(item_settings[item_row,"Parameters"],
                                    getsuffix(item_settings[item_row,"suffix"]),
                                   "=",
                                   item_settings[item_row,"Value"])
                     
                     paras2<-c(paras2,parai2)
                   }
                   items<-paste0(paras2,collapse = "\n")
                   if(item_silent==F) {cat(items)}
                   
                   return(list("themes"=themes,"items"=items))
}




copy_para_csv<-function(output_dir="",
                        para_csvs,
                        from_path=paste0(projectdir,"/03_backups/para_csvs/")){

                          if (!dir.exists(paste0("./",output_dir,"/para_csvs"))) {
                               dir.create(paste0("./",output_dir,"/para_csvs"))
                          }
                          for (para_csv in para_csvs){
                            if (!file.exists(paste0("./",output_dir,"/para_csvs/",para_csv))) {
                                file.copy(from =paste0(from_path,para_csv),
                                          to=   paste0("./",output_dir,"/para_csvs/",para_csv),
                                          copy.mode = F,
                                          copy.date = T)   
                                message(paste0("Copying ",para_csv," to para_csvs folder.\n"))  
                            }else{
                                message(paste0("Found ",para_csv,",use its settings to customise figures.\n"))
                              }
                          }  
  
                       }
