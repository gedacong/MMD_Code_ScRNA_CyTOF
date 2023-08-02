
#' @title imcAsinh
#' @description
#' 处理imc数据的时候使用的转换函数
#'
#' @param value      Vector变量，待转换的数据   
#' @param cofactor   arcsih参数，默认为5

#' @export


imcAsinh<-function (value, cofactor = 5) 
{
  value <- value - 0.05
  loID <- which(value <= 0)
  if (length(loID) > 0) {
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  }
  value <- value/cofactor
  value <- asinh(value)
  return(value)
}