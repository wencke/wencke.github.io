#' 
#' @name GOCluster
#' @title Circular dendrogram.
#' @description GOCluster generates a circular dendrogram of the \code{data} 
#'   clustering using by default euclidean distance and average linkage.The 
#'   inner ring displays the color coded logFC while the outside one encodes the
#'   assigned terms to each gene.
#' @param data A data frame which should be the result of 
#'   \code{\link{circle_dat}} in case the data contains only one logFC column. 
#'   Otherwise \code{data} is a data frame whereas the first column contains the
#'   genes, the second the term and the following columns the logFCs of the 
#'   different contrasts.
#' @param process A character vector of selected processes (ID or term description)
#' @param metric A character vector specifying the distance measure to be used 
#'   (default='euclidean'), see \code{dist}
#' @param clust A character vector specifying the agglomeration method to be 
#'   used (default='average'), see \code{hclust}
#' @param clust.by A character vector specifying if the clustering should be 
#'   done for gene expression pattern or functional categories. By default the 
#'   clustering is done based on the functional categories.
#' @param nlfc If TRUE \code{data} contains multiple logFC columns (default= 
#'   FALSE)
#' @param lfc.col Character vector to define the color scale for the logFC of 
#'   the form c(high, midpoint,low)
#' @param lfc.min Specifies the minimium value of the logFC scale (default = -3)
#' @param lfc.max Specifies the maximum value of the logFC scale (default = 3)
#' @param lfc.space The space between the leafs of the dendrogram and the ring 
#'   for the logFC
#' @param lfc.width The width of the logFC ring
#' @param term.col A character vector specifying the colors of the term bands
#' @param term.space The space between the logFC ring and the term ring
#' @param term.width The width of the term ring
#' @import ggplot2
#' @import ggdendro
#' @import RColorBrewer
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Generating the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Creating the cluster plot
#' GOCluster(circ, EC$process)
#' 
#' #Cluster the data according to gene expression and assigning a different color scale for the logFC
#' GOCluster(circ,EC$process,clust.by='logFC',lfc.col=c('darkgoldenrod1','black','cyan1'))
#' }
#' @export
#' 

GOCluster<-function(data, process, metric, clust, clust.by, nlfc, lfc.col, lfc.min, lfc.max, lfc.space, lfc.width, term.col, term.space, term.width){
  x <- y <- xend <- yend <- width <- space <- logFC <- NULL
  if (missing(metric)) metric<-'euclidean'
  if (missing(clust)) clust<-'average'
  if (missing(clust.by)) clust.by<-'term'
  if (missing(nlfc)) nlfc<-F
  if (missing(lfc.col)) lfc.col<-c('firebrick1','white','dodgerblue')
  if (missing(lfc.min)) lfc.min <- -3
  if (missing(lfc.max)) lfc.max <- 3
  if (missing(lfc.space)) lfc.space<- (-0.5) else lfc.space<-lfc.space*(-1)
  if (missing(lfc.width)) lfc.width<- (-1.6) else lfc.width<-lfc.space-lfc.width-0.1
  if (missing(term.col)) term.col<-brewer.pal(length(process), 'Set3')
  if (missing(term.space)) term.space<- lfc.space+lfc.width else term.space<-term.space*(-1)+lfc.width
  if (missing(term.width)) term.width<- 2*lfc.width+term.space else term.width<-term.width*(-1)+term.space
  
  if (nlfc){
    colnames(data)[1:3]<-c('genes','term','logFC')
    chord<-chord_dat(data[,1:3])
  }else{
    tmp<-data[!duplicated(data$genes),c(5,6)]
    chord<-chord_dat(data,tmp,process)
  }
  if (clust.by=='logFC') distance<-dist(chord[,dim(chord)[2]], method=metric)
  if (clust.by=='term') distance<-dist(chord, method=metric)
  cluster<-hclust(distance, method=clust)
  dendr<-dendro_data(cluster)
  y_range<-range(dendr$segments$y)
  x_pos<-data.frame(x=dendr$label$x, label=as.character(dendr$label$label))
  chord<-as.data.frame(chord)
  chord$label<-as.character(rownames(chord))
  all<-merge(x_pos,chord,by='label')
  all$label<-as.character(all$label)
  if (nlfc){
    lfc_rect<-all[,c(2,dim(all)[2])]
    for (l in 4:dim(data)[2]) lfc_rect<-cbind(lfc_rect, sapply(all$label,function(x) data[match(x,data$genes),l]))
    num<-dim(data)[2]-1
    tmp<-seq(lfc.space,lfc.width,length=num)
    lfc<-data.frame(x=numeric(),width=numeric(),space=numeric(),logFC=numeric())
    for (l in 1:(length(tmp)-1)){
      tmp_df<-data.frame(x=lfc_rect[,1],width=tmp[l+1],space=tmp[l],logFC=lfc_rect[,l+1])
      lfc<-rbind(lfc,tmp_df)
    }
  }else{
    lfc<-all[,c(2,dim(all)[2])]  
    lfc$space<-lfc.space
    lfc$width<-lfc.width
  }
  term<-all[,c(2:(length(process)+2))]
  color<-NULL;termx<-NULL;tspace<-NULL;twidth<-NULL
  for (row in 1:dim(term)[1]){
    idx<-which(term[row,-1]!=0)
    termx<-c(termx,rep(term[row,1],length(idx)))
    color<-c(color,term.col[idx])
    tmp<-seq(term.space,term.width,length=length(idx)+1)
    tspace<-c(tspace,tmp[1:(length(tmp)-1)])
    twidth<-c(twidth,tmp[2:length(tmp)])      
  }
  tmp<-sapply(lfc$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
  logFC<-sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
  lfc$logFC<-logFC
  term_rect<-data.frame(x=termx, width=twidth, space=tspace, col=color)
  legend<-data.frame(x=1:length(process),label=process)

  ggplot()+
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend))+
    geom_rect(data=lfc,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space,fill=logFC))+
    scale_fill_gradient2('logFC',low=lfc.col[3],mid=lfc.col[2],high=lfc.col[1],guide=guide_colorbar(title.position='top',title.hjust=0.5),breaks=c(min(lfc$logFC),max(lfc$logFC)),labels=c(round(min(lfc$logFC)),round(max(lfc$logFC))))+
    geom_rect(data=term_rect,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space),fill=term_rect$col)+
    geom_point(data=legend,aes(x=x,y=0.1,size=factor(label,levels=label),shape=NA))+
    guides(size=guide_legend("GO Terms",ncol=4,byrow=T,override.aes=list(shape=22,fill=term.col,size = 8)))+
    coord_polar()+
    scale_y_reverse()+
    theme(legend.position='bottom',legend.background = element_rect(fill='transparent'),legend.box='horizontal',legend.direction='horizontal')+
    theme_blank  
    
}