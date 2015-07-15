#' circleFun
#' 
#' @name circleFun
#' @param center Center point of the circle
#' @param diameter Defines the diameter of the circle
#' @param npoints Number of points used to draw the circle
#' @import ggplot2
#'   
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
	}

#' get_overlap
#'
#' @name get_overlap
#' @param A Input list 1
#' @param B Input list 2
#' @param C Input list 3
#' @import ggplot2
#' 
get_overlap<-function(A,B,C){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  colnames(C)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(C$logFC)!='numeric'){
    C$logFC<-gsub(",", ".", gsub("\\.", "", C$logFC))
    C$Trend<-sapply(as.numeric(C$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ C$Trend<-sapply(C$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (sum(((A$ID%in%B$ID)==T)==T)==0){
    AB<-data.frame() 
  }else{
    AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    AB<-merge(AB,BA,by="ID")
    rownames(AB)<-AB$ID
    AB<-AB[,-1]   
  }
  if (sum(((A$ID%in%C$ID)==T)==T)==0){
    AC<-data.frame() 
  }else{
    AC<-A[(A$ID%in%C$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    CA<-C[(C$ID%in%A$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    AC<-merge(AC,CA,by="ID")
    rownames(AC)<-AC$ID
    AC<-AC[,-1]
  }
  if (sum(((B$ID%in%C$ID)==T)==T)==0){
    BC<-data.frame() 
  }else{
    BC<-B[(B$ID%in%C$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    CB<-C[(C$ID%in%B$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
    BC<-merge(BC,CB,by="ID")
    rownames(BC)<-BC$ID
    BC<-BC[,-1]
  }
  if (sum(((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T))==0){
    ABC<-data.frame() 
  }else{
    ABC<-A[((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T),which(colnames(A)%in%c('ID','logFC','Trend'))]
    BAC<-B[((B$ID%in%A$ID)==T & (B$ID%in%C$ID)==T),which(colnames(B)%in%c('ID','logFC','Trend'))]
    CAB<-C[((C$ID%in%A$ID)==T & (C$ID%in%B$ID)==T),which(colnames(C)%in%c('ID','logFC','Trend'))]
    ABC<-merge(ABC,BAC,by='ID')
    ABC<-merge(ABC,CAB,by='ID')
    rownames(ABC)<-ABC$ID
    ABC<-ABC[,-1]
  }
  A_only<-A[((A$ID%in%B$ID)==F & (A$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  B_only<-B[((B$ID%in%A$ID)==F & (B$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  C_only<-C[((C$ID%in%A$ID)==F & (C$ID%in%B$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
  rownames(C_only)<-C_only$ID
  C_only<-C_only[,-1]
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  UP<-c(UP,sum(C_only$Trend=='UP'));DOWN<-c(DOWN,sum(C_only$Trend=='DOWN'));Change<-c(Change,sum(C_only$Trend=='Change'))
  if (dim(AB)[1]==0){
    OvAB<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
    OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
    rownames(OvAB)<-rownames(AB)
    AB<-OvAB[order(OvAB$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(AC)[1]==0){
    OvAc<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(AC)[1]) tmp<-c(tmp,ifelse(AC$Trend.x[t]==AC$Trend.y[t],AC$Trend.x[t],'Change'))
    OvAC<-data.frame(logFC_A=AC$logFC.x,logFC_C=AC$logFC.y,Trend=tmp)
    rownames(OvAC)<-rownames(AC)
    AC<-OvAC[order(OvAC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(BC)[1]==0){
    OvBC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(BC)[1]) tmp<-c(tmp,ifelse(BC$Trend.x[t]==BC$Trend.y[t],BC$Trend.x[t],'Change'))
    OvBC<-data.frame(logFC_B=BC$logFC.x,logFC_C=BC$logFC.y,Trend=tmp)
    rownames(OvBC)<-rownames(BC)
    BC<-OvBC[order(OvBC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  if (dim(ABC)[1]==0){
    OvABC<-data.frame()
    UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
  }else{
    tmp<-NULL
    for (t in 1:dim(ABC)[1]) tmp<-c(tmp,ifelse(((ABC$Trend.x[t]==ABC$Trend.y[t]) & (ABC$Trend.x[t]==ABC$Trend[t])),ABC$Trend.x[t],'Change'))
    OvABC<-data.frame(logFC_A=ABC$logFC.x,logFC_B=ABC$logFC.y,logFC_C=ABC$logFC,Trend=tmp)
    rownames(OvABC)<-rownames(ABC)
    ABC<-OvABC[order(OvABC$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  }
  counts<-data.frame(Contrast=c('A_only','B_only','C_only','AB','AC','BC','ABC'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(C_only)[1],dim(AB)[1],dim(AC)[1],dim(BC)[1],dim(ABC)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,C_only=C_only,AB=AB,BC=BC,AC=AC,ABC=ABC)
  return(list(venn_df=counts,table=venn))
}

#' get_overlap2
#'
#' @name get_overlap2
#' @param A Input list 1
#' @param B Input list 2
#' @import ggplot2
#'
get_overlap2<-function(A,B){
  colnames(A)<-c('ID','logFC')
  colnames(B)<-c('ID','logFC')
  UP<-NULL;DOWN<-NULL;Change<-NULL
  if (class(A$logFC)!='numeric'){
    A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
    A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  if (class(B$logFC)!='numeric'){
    B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
    B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN')) 
  }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
  AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
  BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
  A_only<-A[(A$ID%in%B$ID)==F,which(colnames(A)%in%c('ID','logFC','Trend'))]
  B_only<-B[(B$ID%in%A$ID)==F,which(colnames(B)%in%c('ID','logFC','Trend'))]
  AB<-merge(AB,BA,by='ID')
  UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
  UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
  rownames(A_only)<-A_only$ID
  A_only<-A_only[,-1]
  A_only<-A_only[order(A_only$Trend),]
  rownames(B_only)<-B_only$ID
  B_only<-B_only[,-1]
  B_only<-B_only[order(B_only$Trend),]
  tmp<-NULL
  for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
  OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
  rownames(OvAB)<-AB$ID
  AB<-OvAB[order(OvAB$Trend),]
  UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
  counts<-data.frame(Contrast=c('A_only','B_only','AB'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(AB)[1]),UP=UP,DOWN=DOWN,Change=Change)
  venn<-list(A_only=A_only,B_only=B_only,AB=AB)
  return(list(venn_df=counts,table=venn,dim=c(dim(A)[1],dim(B)[1])))
}

#' 
#' @name GOVenn
#' @title Venn diagram of differentially expressed genes.
#' @description The function compares lists of differentially expressed genes 
#'   and illustrates possible relations.Additionally it represents the variety 
#'   of gene expression patterns within the intersection in small pie charts 
#'   with three segements. Clockwise are shown the number of commonly up- 
#'   regulated, commonly down- regulated and contra- regulated genes.
#' @param data1 A data frame consisting of two columns: ID, logFC
#' @param data2 A data frame consisting of two columns: ID, logFC
#' @param data3 A data frame consisting of two columns: ID, logFC
#' @param title The title of the plot
#' @param label A character vector to define the legend keys
#' @param lfc.col A character vector determining the background colors of the 
#'   pie segments representing up- and down- regulated genes
#' @param circle.col A character vector to assign clockwise colors for the 
#'   circles
#' @param plot If TRUE only the venn diagram is plotted. Otherwise the function 
#'   returns a list with two items: the actual plot and a list containing the 
#'   overlap entries (default= TRUE)
#' @details The \code{plot} argument can be used to adjust the amount of 
#'   information that is returned by calling the function. If you are only 
#'   interested in the actual plot of the venn diagram, \code{plot} should be 
#'   set to TRUE. Sometimes you also want to know the elements of the 
#'   intersections. In this case \code{plot} should be set to FALSE and the 
#'   function call will return a list of two items. The first item, that can be 
#'   accessed by $plot, contains the plotting information. Additionally, a list
#'   ($table) will be returned containing the elements of the various overlaps. 
#' @import ggplot2
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Generating the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Selecting terms of interest
#' l1<-subset(circ,term=='heart development',c(genes,logFC))
#' l2<-subset(circ,term=='plasma membrane',c(genes,logFC))
#' l3<-subset(circ,term=='tissue morphogenesis',c(genes,logFC))
#' 
#' GOVenn(l1,l2,l3, label=c('heart development','plasma membrane','tissue morphogenesis'))
#' }
#' @export

GOVenn<-function(data1, data2, data3, title, label, lfc.col, circle.col, plot=T){
  id <- NULL
  if (missing(label)) label<-c('List1','List2','List3')
  if (missing(lfc.col)) lfc.col<-c('firebrick1','gold','cornflowerblue')
  if (missing(circle.col)) circle.col<-c('brown1','chartreuse3','cornflowerblue')
  if (missing(title)) title<-''
  if (missing(data3)==F) {
    three<-T
    overlap<-get_overlap(data1,data2,data3)
    venn_df<-overlap$venn_df
    table<-overlap$table
  }else{
    three<-F
    overlap<-get_overlap2(data1,data2)
    venn_df<-overlap$venn_df
    table<-overlap$table
  }
  
	### calc Venn ###
  if (three){
    center<-data.frame(x=c(0.4311,0.4308,0.6380),y=c(0.6197,0.3801,0.5001),diameter=c(0.4483,0.4483,0.4483))
    outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
	  for (var in 1:3){
		  dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
		  outerCircle<-rbind(outerCircle,dat)
	  }
	  outerCircle$id<-rep(c(label[1],label[2],label[3]),each=100)
	  outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2],label[3]))
  }else{
    center<-data.frame(x=c(0.33,0.6699),y=c(0.5,0.5),diameter=c(0.6180,0.6180))
    outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
    for (var in 1:2){
      dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
      outerCircle<-rbind(outerCircle,dat)
    }
    outerCircle$id<-rep(c(label[1],label[2]),each=100)
    outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2]))
  }

	### calc single pies ### 
  if (three){
    Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
	  dat <- circleFun(c(center$x[1],max(subset(outerCircle,id==label[1])$y)-0.05),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  dat <- circleFun(c(center$x[2],min(subset(outerCircle,id==label[2])$y)+0.05),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  dat <- circleFun(c(max(subset(outerCircle,id==label[3])$x)-0.05,center$y[3]),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  Pie$id<-rep(1:3,each=100)
	  UP<-Pie[c(1:50,100:150,200:250),]
	  Down<-Pie[c(50:100,150:200,250:300),]
  }else{
    Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
    dat <- circleFun(c(min(subset(outerCircle,id==label[1])$x)+0.05,center$y[1]),0.1,npoints = 100)
    Pie<-rbind(Pie,dat)
    dat <- circleFun(c(max(subset(outerCircle,id==label[2])$x)-0.05,center$y[2]),0.1,npoints = 100)
    Pie<-rbind(Pie,dat)
    Pie$id<-rep(1:2,each=100)
    UP<-Pie[c(1:50,100:150),]
    Down<-Pie[c(50:100,150:200),]
  }

	### calc single pie text ###
	if (three){
    x<-c();y<-c()
	  for (i in unique(Pie$id)){
		  x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
		  y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
		  y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
	  }
	  pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2],venn_df$UP[3],venn_df$DOWN[3]))
  }else{
    x<-c();y<-c()
    for (i in unique(Pie$id)){
      x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
      y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
      y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
    }
    pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2]))
  }
  
	### calc overlap pies ### 
  if (three){
    smc<-data.frame(x=c(0.6,0.59,0.31,0.5),y=c(0.66,0.34,0.5,0.5))
    PieOv<-data.frame(x=numeric(),y=numeric())
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[1],smc$y[1]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[2],smc$y[2]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[3],smc$y[3]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[4],smc$y[4]),0.06,npoints = 100))
	  PieOv$id<-rep(1:4,each=100)
	  smc$id<-1:4
	  UPOv<-rbind(smc[1,],PieOv[1:33,],smc[1,],smc[2,],PieOv[100:133,],smc[2,],smc[3,],PieOv[200:233,],smc[3,],smc[4,],PieOv[300:333,],smc[4,])
	  Change<-rbind(smc[1,],PieOv[33:66,],smc[1,],smc[2,],PieOv[133:166,],smc[2,],smc[3,],PieOv[233:266,],smc[3,],smc[4,],PieOv[333:366,],smc[4,])
	  DownOv<-rbind(smc[1,],PieOv[66:100,],smc[1,],smc[2,],PieOv[166:200,],smc[2,],smc[3,],PieOv[266:300,],smc[3,],smc[4,],PieOv[366:400,],smc[4,])
  }else{
    PieOv<-data.frame(x=numeric(),y=numeric(),id=numeric())
    PieOv<-rbind(PieOv,circleFun(c(0.5,0.5),0.08,npoints = 100))
    PieOv$id<-rep(1,100)
    center<-data.frame(x=0.5, y=0.5, id=1)
    UPOv<-rbind(center[1,],PieOv[1:33,])
    Change<-rbind(center[1,],PieOv[33:66,])
    DownOv<-rbind(center[1,],PieOv[66:100,])
  }
  
  ### calc overlap pie text ###
  if (three){
    x<-c();y<-c()
	  for (i in unique(PieOv$id)){
		  x<-c(x,subset(UPOv,id==i)$x[1]+0.0115,subset(DownOv,id==i)$x[1]-0.018,subset(Change,id==i)$x[1]+0.01)
		  y<-c(y,subset(UPOv,id==i)$y[1]+0.01,subset(DownOv,id==i)$y[1],subset(Change,id==i)$y[1]-0.013)
	  }
	  small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[5],venn_df$Change[5],venn_df$DOWN[5],venn_df$UP[6],venn_df$Change[6],venn_df$DOWN[6],venn_df$UP[4],venn_df$Change[4],venn_df$DOWN[4],venn_df$UP[7],venn_df$Change[7],venn_df$DOWN[7]))
  }else{
    x<-c(subset(UPOv,id==1)$x[1]+0.015,subset(DownOv,id==1)$x[1]-0.018,subset(Change,id==1)$x[1]+0.01)
    y<-c(subset(UPOv,id==1)$y[1]+0.015,subset(DownOv,id==1)$y[1],subset(Change,id==1)$y[1]-0.013)
    small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[3],venn_df$Change[3],venn_df$DOWN[3]))
  }
  
	g<- ggplot()+
    geom_polygon(data=outerCircle, aes(x,y, group=id, fill=id) ,alpha=0.5,color='black')+
	  scale_fill_manual(values=circle.col)+
    guides(fill=guide_legend(title=''))+
	  geom_polygon(data=UP, aes(x,y,group=id),fill=lfc.col[1],color='white')+
	  geom_polygon(data=Down, aes(x,y,group=id),fill=lfc.col[3],color='white')+
	  geom_text(data=pieText, aes(x=x,y=y,label=label),size=5)+
    geom_polygon(data=UPOv, aes(x,y,group=id),fill=lfc.col[1],color='white')+
	  geom_polygon(data=DownOv, aes(x,y,group=id),fill=lfc.col[3],color='white')+
	  geom_polygon(data=Change, aes(x,y,group=id),fill=lfc.col[2],color='white')+
	  geom_text(data=small.pieT,aes(x=x,y=y,label=label),size=4)+
    theme_blank+
    labs(title=title)
  
  if (plot) return(g) else return(list(plot=g,table=table))
}


	
