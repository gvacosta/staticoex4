#' Table preprocessing: Centering and Normalizing
#'
#' @description
#' Scaling of two dataframes: environmental and species variables. There are three options for preprocessing of tables regarding centering and normalizing.
#'
#' @usage Prep.Tables(x, y, condition, preprocessing)
#'
#' @param x A dataframe with k matrices juxtaposed in vertical order containing the data of the P variables (environmental information).
#' @param y A dataframe with k matrices juxtaposed in vertical order containing the data of the Q variables (species information).
#' @param condition Factor, the k condition descriptions that will define the third path for the analysis.
#' @param preprocessing A string indicating the standardization method to be used: "B.Total" performs total standardization and "B.Partial" partial standardization of Bouroche (1975) (library ade4) and "A.Norma" standardization of Abdi et al. (2012). See details.
#'
#' @details
#' Table preprocessing function that considers a third way or categorical variable that divides data into k tables. That is, the k pairs of individual tables must be appended both in the matrix of the environmental variables and in the species data and, in both cases, have the same number of observations, since an index table is built and runs through the same positions of the observations for both matrices. There are three options for matrix scaling: centering and normalizing to remove the effects of different means and scales.
#' •	“B.Partial” center and normalize the columns for each k section, ensuring that the initial k matrices have zero mean and unit variance.
#' •	“B.Total” center the columns for each k section, then normalize for the total column, ensuring that the initial k-matrices have a zero mean, and the variance of the total column is one.
#' •	“A.Norma” This method operates similarly to dividing each element of the total matrix by the square root of the total sum of squared elements, ensuring that the total variance of each k-matrix is equal to one and that all K matrices participate equally in the analysis.
#' It is worth mentioning that X is normalized, and Y is only centered, due to X is taken as variables or measurements with possible different scales and Y as frequency data (species frequencies).
#'
#' @return Returns a list-type object with three elements:
#' \itemize{
#' \item *X:* It is the normalized matrix corresponding to “x” (environmental matrix).
#' \item *Y:* It is the **centered** matrix corresponding to “y” (species matrix).
#' \item *Index:* A data.frame with the name of k-tables, the number of observations and the index of start and finish for each table.
#' }
#'
#' @author Mariela Gonzalez Narvaez, Ph.D., María José Fernández Gomez, Susana Mendes, Purificación Galindo-Villardón, Ph.D., Gema Zambrano Zambrano, Ysaí Ronquillo Mora, George Acosta Chong
#'
#' @references http://hdl.handle.net/10366/149381
#'
#' @seealso [HJSTATICO()], ade4::statico, MixSTATICO()
#'
#' @importFrom methods is
#'
#' @examples
#' Prep.Tables(x=env_data, y=spe_data, condition=etiq_data$space, preprocessing = "A.Norma")
#'
#' @export
Prep.Tables=function(x,y,condition,preprocessing) {
  if(exists("preprocessing")==FALSE){
    preprocessing="A.Norma"
  }
  if (is(condition,"factor")){
    condition<-as.character(condition)
  }
  dimenx=dim(x); dimeny=dim(y)
  k=length(unique(condition))
  #Index matrix environment and species
  index=matrix(nrow=k,ncol=4)
  index=as.data.frame(index)
  colnames(index)=c("Tables","Start","Finish","n")
  index[1,1]=condition[1]
  index[1,2]=1
  count=1##;i=1 #no es necesario definir i
  s=1; num=1
  for (i in 2:length(condition)) {
    if(condition[i]!=index[s,1]) {
      index[s,3]=count
      index[s,4]=num
      s=s+1
      num=1
      index[s,1]=condition[i]
      count=count+1
      index[s,2]=count}
    else {
      count=count+1
      num=num+1
      index[s,3]=count
      index[s,4]=num
    }
  }
  index
  #standarization by column - whole table
  mean=apply(x,2,mean)
  sd=apply(x,2,sd)
  sd= sd*sqrt((dimenx[1]-1)/dimenx[1])
  z=matrix(nrow=dimenx[1],ncol=dimenx[2])

  for (j in 1:dimenx[2]) {
    z[,j]=(x[,j]-mean[j])/sd[j]
  }
  if (preprocessing=="A.Norma") {
    #preprocesamiento de las tablas,centrar y normalizar cada tabla por columnas
    for(i in 1:k){
      temp=matrix(nrow=index[i,4],ncol=dimenx[2])
      temp=z[index[i,2]:index[i,3],]
      temp=as.matrix(temp)
      temp2=temp%*%t(temp)
      sol=1/sqrt(norm(temp2,type=c("O")))#norma
      if (i==1) {
        X=z[index[i,2]:index[i,3],]*sol}
      else{
        temp3=z[index[i,2]:index[i,3],]*sol
        X=rbind(X,temp3)
      }
    }
  }
  if (preprocessing=="B.Partial") {
    avem=matrix(nrow=dimenx[1],ncol=dimenx[2])
    stdm=matrix(nrow=dimenx[1],ncol=dimenx[2])
    for (i in 1:k) {
      temp=z[index[i,2]:index[i,3],]
      tempav=apply(temp,2,mean)
      tempsd=apply(temp,2,sd)
      tempsd=tempsd*sqrt((index[i,4]-1)/index[i,4])
      ta=matrix(nrow=index[i,4],ncol=dimenx[2])
      tsd=matrix(nrow=index[i,4],ncol=dimenx[2])
      for (ni in 1:index[i,4]) {
        ta[ni,]=tempav
        tsd[ni,]=tempsd
      }
      if (i==1){
        avem=ta
        stdm=tsd
      }else{
        avem=rbind(avem,ta)
        stdm=rbind(stdm,tsd)
      }
    }
    X=(z-avem)/stdm
  }
  if (preprocessing=="B.Total") {
    avem=matrix(nrow=dimenx[1],ncol=dimenx[2])
    for (i in 1:k) {
      temp=z[index[i,2]:index[i,3],]
      tempav=apply(temp,2,mean)
      ta=matrix(nrow=index[i,4],ncol=dimenx[2])
      for (ni in 1:index[i,4]) {
        ta[ni,]=tempav }
      if (i==1){ avem=ta}
      else { avem=rbind(avem,ta)}
    }
    r=z-avem
    temp=apply(r,2,sd)
    temp=temp*sqrt((dimenx[1]-1)/dimenx[1])
    stdm=matrix(nrow=dimenx[1],ncol=dimenx[2])
    for (i in 1:dimenx[1]) {
      stdm[i,]=temp }
    X=r/stdm
  }
  #Species Matrix-average matrix
  avematrixs=matrix(nrow=dimeny[1],ncol=dimeny[2])
  mean_spe=apply(y,2,mean)
  for (i in 1:dimeny[2]) {
    avematrixs[,i]=mean_spe[i]
  }
  Y=y-avematrixs
  colnames(X)=colnames(x)
  Lst=list(X=X,Y=Y,index=index)
}
