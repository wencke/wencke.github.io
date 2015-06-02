#' Transcriptomic information of endothelial cells.
#' 
#' The data set contains the transcriptomic information of endothelial cells
#' from two steady state tissues (brain and heart). More detailed information
#' can be found in the paper by Nolan et al. 2013. The data was normalized and a
#' statistical analysis was performed to determine differentially expressed
#' genes. DAVID functional annotation tool was used to perform a gene-
#' annotation enrichment analysis of the set of differentially expressed genes
#' (adjusted p-value < 0.05).
#' 
#' @docType data
#' @keywords datasets
#' @name EC
#' @usage data(EC)
#' @format A list containing 5 items 
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47067}
"EC"

#' theme blank
#'
#' @name theme_blank
#' @import ggplot2
#' 

theme_blank <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())

#' draw table
#'
#' @name draw_table
#' @param data data frame; GO ID and GO term
#' @param size numeric; defines font size
#' @param col character; defines color of table header
#' @import ggplot2
#' @import grid
#' @import gridExtra  
#' 

draw_table <- function(data, size = 4, col = ''){
  x <- y <- idx <- NULL
	colnames(data) <- tolower(colnames(data))
  width <- size * 0.0125
	nrow <- dim(data)[1]
	TLen <- max(nchar(data$term))
	IDLen <- max(nchar(data$id))
	xmin <- 0 - (IDLen / 2 + 1) * width; x1 <- 0 + (IDLen / 2 + 1) * width; xmax <- x1 + (TLen + 2) * width
	ymin <- 0 - width; ymax <- nrow * (2 * width); yheader <- (nrow + 1) * (2 * width)
	table <- data.frame(x = c(xmin, x1, x1, xmin, xmin, x1, xmax, xmax, x1, x1), y = rep(c(ymin, ymin, ymax, ymax, ymin), 2), idx = rep(c(1, 2), each = 5))
	text <- data.frame(x = rep(c(0, (x1 + xmax) / 2), each = nrow + 1), y = rep(seq(yheader - width, 0, length = nrow + 1), 2),label = c('GO ID', data$id, 'Term', data$term))
	header <- data.frame(x = c(xmin, x1, x1, xmin, xmin, x1, xmax, xmax, x1, x1), y = rep(c(ymax, ymax, yheader, yheader, ymax), 2), idx = rep(c(1, 2), each = 5))
	if (length(col) == 1){
	  text.col <- rep('black', 2 * nrow + 2)
	}else{
	  text.col <- rep(c('black', rep(col[1], sum(data$category == 'BP')), rep(col[2], sum(data$category == 'CC')), rep(col[3], sum(data$category == 'MF'))), 2)
	}
  g<-	ggplot()+
		geom_polygon(data = table, aes(x, y, group = idx), fill = 'transparent', color = 'black')+
    geom_polygon(data = header, aes(x, y, group = idx), fill = 'gray70', color = 'black')+
    geom_text(data = text, aes(x = x, y = y, label = label), color = text.col, size = rep(4, 2 * nrow + 2))+
    ylim(ymin - 1.5, yheader + 1.5)+
    xlim(xmin - width, xmax + width)+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, -1.1), 'cm'))+
		theme_blank
	return(g)
}

#' bezier
#'
#' @name bezier
#' @param data data frame x.start,x.end,y.start,y.end
#' @param process.col Color of the processes
#' 

bezier <- function(data, process.col){
	x <- c()
	y <- c()
	Id <- c()
	sequ <- seq(0, 1, by = 0.01)
	N <- dim(data)[1]
	sN <- seq(1, N, by = 2)
	if (process.col[1] == '') col_rain <- rainbow(N) else col_rain <- process.col
	for (n in sN){
		xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
		for (t in sequ){
			xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
			xval <- c(xval, xva)
			xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
			xval2 <- c(xval2, xva2)
			yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]	
			yval <- c(yval, yva)
			yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
			yval2 <- c(yval2, yva2)			
		}
		x <- c(x, xval, rev(xval2))
		y <- c(y, yval, rev(yval2))
		Id <- c(Id, rep(n, 2 * length(sequ)))
	}
	df <- data.frame(lx = x, ly = y, ID = Id)
	return(df)
}

#' 
#' @name circle_dat
#' @title Creates a plotting object.
#' @description The function takes the results from a functional analysis (for 
#'   example DAVID) and combines it with a list of selected genes and their 
#'   logFC. The resulting data frame can be used as an input for various ploting
#'   functions.
#' @param terms A data frame with columns for 'category', 'ID', 'term', adjusted
#'   p-value ('adj_pval') and 'genes'
#' @param genes A data frame with columns for 'ID', 'logFC'
#' @details Since most of the gene- annotation enrichment analysis are based on 
#'   the gene ontology database the package was build with this structure in 
#'   mind, but is not restricted to it. Gene ontology is structured as an 
#'   acyclic graph and it provides terms covering different areas. These terms 
#'   are grouped into three independent \code{categories}: BP (biological 
#'   process), CC (cellular component) or MF (molecular function).
#'   
#'   The "ID" and "term" columns of the \code{terms} data frame refer to the ID 
#'   and term description, whereas the ID is optional.
#'   
#'   The "ID" column of the \code{genes} data frame can contain any unique 
#'   identifier. Nevertheless, the identifier has to be the same as in "genes" 
#'   from \code{terms}.
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Building the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' }
#' @export

circle_dat <- function(terms, genes){

	colnames(terms) <- tolower(colnames(terms))
	terms$genes <- toupper(terms$genes)
	genes$ID <- toupper(genes$ID)
	tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
	count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
	logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
	if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
	}
	s <- 1; zsc <- c()
	for (c in 1:length(count)){
		value <- 0
		e <- s + count[c] - 1
		value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
		zsc <- c(zsc, sum(value) / sqrt(count[c]))
		s <- e + 1
	}
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                   count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                   zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
			   count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
			   zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
	return(df)
}

#' 
#' @name chord_dat
#' @title Creates a binary matrix.
#' @description The function creates a matrix which represents the binary 
#'   relation (1= is related to, 0= is not related to) between selected genes 
#'   (row) and processes (column). The resulting matrix can be visualized with 
#'   the \code{\link{GOChord}} function.
#' @param data A data frame with coloumns for GO ID|term and genes.
#' @param genes A character vector of selected genes OR data frame with coloumns
#'   for gene ID and logFC.
#' @param limit A vector with two cutoff values (default= c(0,0)). The first 
#'   value defines the minimum number of terms a gene has to be assigned to. The
#'   second the minimum number of genes assigned to a selected term.
#' @details If more than one logFC value for each gene is at disposal, only one 
#'   should be used to create the binary matrix. The other values have to be 
#'   added manually later. The parameter \code{limit} can be used to reduce the 
#'   dimension of the calculated matrix. This might be useful to represent the 
#'   data more clearly with \code{GOChord} later on. The first value of the 
#'   vector defines the threshold for the minimum number of terms a gene has to 
#'   be assigned to in order to be represented in the plot. Most of the time it 
#'   is more meaningful to represent genes with various functions. A value of 3 
#'   excludes all genes with less than three term assignments. Whereas the 
#'   second value of the parameter restricts the number of terms according to 
#'   the number of assigned genes. All terms with a count smaller or equal to
#'   the threshold are excluded.
#' @param process A character vector of selected processes
#' @return A binary matrix
#' @seealso \code{\link{GOChord}}
#' @examples
#' \dontrun{
#' # Load the included dataset
#' data(EC)
#' 
#' # Building the circ object
#' circ <- circular_dat(EC$david, EC$genelist)
#' 
#' # Building the binary matrix
#' chord <- chord_dat(circ, EC$genes, EC$process)
#' 
#' # Excluding genes which are assigned only to a single term
#' chord <- chord_dat(circ, EC$genes, EC$process, limit = c(1,0))
#' 
#' # Excluding terms with a count smaller than 5
#' chord <- chord_dat(circ, EC$genes, EC$process, limit = c(0,5))
#' 
#' }
#' @export

chord_dat <- function(data, genes, process, limit){
  id <- term <- BPprocess <- NULL
  colnames(data) <- tolower(colnames(data))
  if (missing(limit)) limit <- c(0, 0)
  if (missing(genes)){
    if (is.null(data$logFC)){
      genes <- unique(data$genes)
    }else{
      sub <- data[!duplicated(data$genes), ]
      genes <- data.frame(genes = sub$genes, logFC = sub$logFC)
    }
  }else{
    if(is.vector(genes)){
      genes <- as.character(genes) 
    }else{
      if(class(genes[, 2]) != 'numeric') logFC <- as.numeric(levels(genes[, 2]))[genes[, 2]] else logFC <- genes[,2]
      genes <- as.character(genes[, 1])      
    }
  }
  if (missing(process)){
    process <- unique(data$term)
  }else{
    if(class(process) != 'character') process <- as.character(process)
  }
  if (strsplit(process[1],':')[[1]][1] == 'GO'){
    subData <- subset(data, id%in%process)
    colnames(subData)[which(colnames(subData) == 'id')] <- 'BPprocess'
  }else{
    subData <- subset(data, term%in%process)
    colnames(subData)[which(colnames(subData) == 'term')] <- 'BPprocess'
  }
	
	mat <- matrix(0, ncol = length(process), nrow = length(genes))
	if (length(unique(subData$BPprocess)) < dim(subData)[1]){
		for (p in 1:length(process)){
			sub2 <- subset(subData, BPprocess == process[p])
			for (g in 1:length(genes))mat[g, p] <- ifelse(genes[g]%in%sub2$genes, 1, 0)
		}
		rownames(mat) <- genes
		lfc <- sapply(genes, function(x) subData[match(x, subData$genes), which(colnames(subData) == 'logfc')])
		mat <- cbind(mat, lfc)
		colnames(mat) <- c(process, 'logFC')
	}else{
		row <- 0
		for (g in genes){
			row <- row + 1
			for(p in 1:dim(subData)[1]){
				d <- strsplit(as.vector(subData$genes[p]), ',')
				vec <- c()
				for (sub in 1:length(d[[1]])) vec <- c(vec, ifelse(g == (d[[1]][sub]), 1, 0)) 
				mat[row, p] <- ifelse(sum(vec) == 0, 0, 1)
			}
		}
		if (!is.vector(genes)) mat <- cbind(mat, logFC)
		rownames(mat) <- genes
		if(dim(mat)[2] != length(subData$BPprocess)) colnames(mat) <- c(as.character(subData$BPprocess), 'logFC') else colnames(mat) <- subData$BPprocess		
	}
	idx <- which(colnames(mat) == 'logFC')
  mat <- mat[is.na(mat[, idx]) == F, ]
  if ((sum(rowSums(mat[, -dim(mat)[2]]) < limit[1])) > 0) mat <- mat[-which(rowSums(mat[, -dim(mat)[2]]) < limit[1]), ]
  if ((sum(colSums(mat) < limit[2])) > 0) mat <- mat[, -which(colSums(mat) < limit[2])]
  if ((sum(rowSums(mat[, -dim(mat)[2]]) < limit[1])) > 0) mat <- mat[-which(rowSums(mat[, -dim(mat)[2]]) < limit[1]), ]
  
  return(mat)
}

#' 
#' @name GOBubble
#' @title Bubble plot.
#' @description The function creates a bubble plot of the input \code{data}. The
#'   input \code{data} can be created with the help of the 
#'   \code{\link{circle_dat}} function.
#' @param data A data frame with coloumns for category, GO ID, term, adjusted 
#'   p-value, z-score, count(num of genes)
#' @param display A character vector. Indicates whether it should be a single 
#'   plot ('single') or a facet plot with panels for each category 
#'   (default='single')
#' @param title The title (on top) of the plot
#' @param color A character vector which defines the color of the bubbles for 
#'   each category
#' @param labels Sets a threshold for the displayed labels. The threshold refers
#'   to the -log(adjusted p-value) (default=5)
#' @param ID If TRUE then labels are IDs else terms
#' @param table.legend Defines whether a table of GO ID and GO term should be 
#'   displayed on the right side of the plot or not (default=TRUE)
#' @param table.col If TRUE then the table entries are colored according to 
#'   their category, if FALSE then entries are black
#' @param bg If TRUE within the plot a grey background will be drawn
#' @details The x- axis of the plot represents the z-score. The negative
#'   logarithm of the adjusted p-value (corresponding to the significance of the
#'   term) is displayed on the y-axis. The area of the plotted circles is 
#'   proportional to the number of genes assigned to the term. Each circle is 
#'   colored according to its category and labeled alternatively with the ID or 
#'   term name.
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import ggvis
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Building the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Creating the bubble plot coloring the table entries according to the category
#' GOBubble(circ, table.col=T)
#' 
#' #Creating the bubble plot displaying the term instead of the ID and without the table
#' GOBubble(circ,ID=F,table.legend=F)
#' 
#' #Faceting the plot
#' GOBubble(circ, display='multiple')
#' }
#' @export

GOBubble <- function(data, display, title, color, labels, ID = T, table.legend = T, table.col = T, static = T){
  zscore <- adj_pval <- category <- count <- id <- term <- NULL
  if (missing(display)) display <- 'single'
  if (missing(title)) title <- ''
  if (missing(color)) cols <- c("chartreuse4", "brown2", "cornflowerblue") else cols <- color
  if (missing(labels)) labels <- 5
  
  if(static){
    colnames(data) <- tolower(colnames(data))
	  if(!'count'%in%colnames(data)){
	    rang <- c(5, 5)
	    data$count <- rep(1, dim(data)[1])
	  }else {rang <- c(1, 30)}
	  data$adj_pval <- -log(data$adj_pval, 10)
	  sub <- data[!duplicated(data$term), ]
	  g <- ggplot(sub, aes(zscore, adj_pval))+
      labs(title = title, x = 'z-score', y = '-log (adj p-value)')+
      geom_point(aes(col = category, size = count), alpha = 1 / 2)+
      geom_hline(yintercept = 1.3, col = 'orange')+
      scale_size(range = rang, guide = 'none')
	  if (!is.character(labels)) sub2 <- subset(sub, subset = sub$adj_pval >= labels) else sub2 <- subset(sub, sub$id%in%labels | sub$term%in%labels)
	  if (display == 'single'){
		  g <- g + scale_colour_manual('Category', values = cols, labels = c('Biological Process', 'Cellular Component', 'Molecular Function'))+
		    theme(legend.position = 'bottom')+
		    annotate ("text", x = min(sub$zscore), y = 1.5, label = "threshold", colour = "orange", size = 3)
		  if (ID) g <- g+ geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5) else g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 4)
      if (table.legend){
        if (table.col) table <- draw_table(sub2, size = 4, col = cols) else table <- draw_table(sub2, size = 4, col = 'black')
        g <- g + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), 'cm'), axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.background = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major = element_line(color = 'grey80'), plot.background = element_blank()) 
        par(mar = c(0.1, 0.1, 0.1, 0.1))
			  grid.arrange(g, table, heights = c(1, 1), widths = c(1.8, 1), ncol = 2)
		  }else{
        g + theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.background = element_blank(),
                  panel.grid.minor = element_blank(), panel.grid.major = element_line(color = 'grey80'), plot.background = element_blank(),
                  axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'))
      }
	  }else{
		  g <- g + facet_grid(.~category, space = 'free_x', scales = 'free_x') + scale_colour_manual(values = cols, guide ='none')
		  if (ID) {
          g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5) + 
            theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.border = element_rect(fill = 'transparent', color = 'grey80'),
                  panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank()) 
      }else{
          g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 5) + 
            theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.border = element_rect(fill = 'transparent', color = 'grey80'),
                  panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank())
      }
	  }
 }else{
   colnames(data) <- tolower(colnames(data))
   data$adj_pval <- -log(data$adj_pval, 10)
   sub <- data[!duplicated(data$term), ]
   if(!'count'%in%colnames(data)){
     rang <- c(400, 400)
     data$count <- rep(1, dim(data)[1])
   }else {rang <- c(50, 800)}
   texts <- function(x){
     if(is.null(x)) return(NULL)
     paste0("GO Term: ", x$term)
   }
   sub %>% ggvis(~zscore, ~adj_pval, fill = ~category, size = ~count, key := ~term, opacity := 0.5) %>% 
     layer_points() %>%
     scale_numeric('size', range = rang) %>%
     hide_legend('fill') %>%
     add_tooltip(texts, 'hover')
 }
}

#' 
#' @name GOBar
#' @title Z-score colored barplot.
#' @description Z-score colored barplot of terms ordered alternatively by 
#'   z-score or the negative logarithm of the adjusted p-value
#' @param data A data frame containing at least the term ID and/or term, the 
#'   adjusted p-value and the z-score. A possible input can be generated with 
#'   the \code{circle_dat} function
#' @param display A character vector indicating whether a single plot ('single')
#'   or a facet plot with panels for each category should be drawn 
#'   (default='single')
#' @param order.by.zscore Defines the order of the bars. If TRUE the bars are 
#'   ordered according to the z-scores of the processes. Otherwise the bars are 
#'   ordered by the negative logarithm of the adjusted p-value
#' @param title The title of the plot
#' @param zsc.col Character vector to define the color scale for the z-score of 
#'   the form c(high, midpoint,low)
#' @details If \code{display} is used to facet the plot the width of the panels
#'   will be proportional to the length of the x scale.
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import ggvis
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Building the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Creating the bar plot
#' GOBar(circ)
#' 
#' #Faceting the plot
#' GOBar(circ, display='multiple')
#' }
#' @export

GOBar <- function(data, display, order.by.zscore = T, title, zsc.col, static = T){
  id <- adj_pval <- zscore <- NULL
  if (missing(display)) display <- 'single'
  if (missing(title)) title <- ''
  if (missing(zsc.col)) zsc.col <- c('red', 'white', 'blue')
  colnames(data) <- tolower(colnames(data))
	data$adj_pval <- -log(data$adj_pval, 10)
	sub <- data[!duplicated(data$term), ]
  
  if (static){ 
	  if (order.by.zscore == T) {
		  sub <- sub[order(sub$zscore, decreasing = T), ]
		  leg <- theme(legend.position = 'bottom')
		  g <-  ggplot(sub, aes(x = factor(id, levels = reorder(id, adj_pval)), y = adj_pval, fill = zscore)) +
		    geom_bar(stat = 'identity', color = 'black') +
		    scale_fill_gradient2('z-score', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), 
                             breaks = c(min(sub$zscore), max(sub$zscore)), labels = c('decreasing', 'increasing')) +
		    labs(title = title, x = '', y = '-log (adj p-value)') +
		    leg
	  }else{
		  sub <- sub[order(sub$adj_pval, decreasing = T), ]
		  leg <- theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.995), legend.background = element_rect(fill = 'transparent'),
                   legend.box = 'vertical', legend.direction = 'horizontal')
		  g <-  ggplot(sub, aes( x = factor(id, levels = reorder(id, adj_pval)), y = zscore, fill = adj_pval)) +
		    geom_bar(stat = 'identity', color = 'black') +
		    scale_fill_gradient2('Significance', guide = guide_colorbar(title.position = "top", title.hjust = 0.5), breaks = c(min(sub$adj_pval), max(sub$adj_pval)), labels = c('low', 'high')) +
		    labs(title = title, x = '', y = 'z-score') +
		    leg
	  }
	  if (display == 'single'){
      g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'),
               panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), plot.background = element_blank())        
    }else{
      g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'),
               panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), plot.background = element_blank())+
               facet_grid(.~category, space = 'free_x', scales = 'free_x')
    }
  }else{
    texts <- function(x){
      if(is.null(x)) return(NULL)
      paste0("GO Term: ", x$term)
    }
    sub %>% 
      ggvis(~id, ~adj_pval, key:= ~term) %>% 
      layer_bars() %>%
      add_tooltip(texts, 'hover')
  }
}

#' 
#' @name GOCircle
#' @title Circular visualization of the results of a functional analysis.
#' @description The circular plot combines gene expression and gene- annotation 
#'   enrichment data. A subset of terms is displayed like the \code{GOBar} plot 
#'   in combination with a scatterplot of the gene expression data. The whole 
#'   plot is drawn on a specific coordinate system to achieve the circular 
#'   layout.The segments are labeled with the term ID.
#' @param data A special data frame which should be the result of 
#'   \code{circle_dat}
#' @param title The title of the plot
#' @param nsub A numeric or character vector. If it's numeric then the number 
#'   defines how many processes are displayed (starting from the first row of 
#'   \code{data}). If it's a character string of processes then these processes 
#'   are displayed
#' @param rad1 The radius of the inner circle (default=2)
#' @param rad2 The radius of the outer circle (default=3)
#' @param table.legend Shall a table be displayd or not? (default=TRUE)
#' @param zsc.col Character vector to define the color scale for the z-score of 
#'   the form c(high, midpoint,low)
#' @param lfc.col A character vector specifying the color for up- and 
#'   down-regulated genes
#' @param label.size Size of the segment labels (default=5)
#' @param label.fontface Font style of the segment labels (default='bold')
#' @details The outer circle shows a scatter plot for each term of the logFC of 
#'   the assigned genes. The colors can be changed with the argument 
#'   \code{lfc.col}.
#'   
#'   The \code{nsub} argument needs a bit more explanation to be used wisely. First of 
#'   all, it can be a numeric or a character vector. If it is a character vector
#'   then it contains the IDs or term descriptions of the displayed processes.If
#'   \code{nsub} is a numeric vector then the number defines how many terms are 
#'   displayed. It starts with the first row of the input data frame.
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @seealso \code{\link{circle_dat}}, \code{\link{GOBar}}
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Building the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Creating the circular plot
#' GOCircle(circ)
#' 
#' #Creating the circular plot with a different color scale
#' GOCircle(circ,zsc.col=c('yellow','black','cyan'))
#' 
#' #Creating the circular plot with different font style
#' GOCircle(circ,label.size=5,label.fontface='italic')
#' }
#' @export

GOCircle <- function(data, title, nsub, rad1, rad2, table.legend = T, zsc.col, lfc.col, label.size, label.fontface){
  xmax <- y1<- zscore <- y2 <- ID <- logx <- logy2 <- logy <- logFC <- NULL
  if (missing(title)) title <- ''
  if (missing(nsub)) if (dim(data)[1] > 10) nsub <- 10 else nsub <- dim(data)[1]
  if (missing(rad1)) rad1 <- 2
  if (missing(rad2)) rad2 <- 3
  if (missing(zsc.col)) zsc.col <- c('red', 'white', 'blue')
  if (missing(lfc.col)) lfc.col <- c('cornflowerblue', 'firebrick1') else lfc.col <- rev(lfc.col)
  if (missing(label.size)) label.size = 5
  if (missing(label.fontface)) label.fontface = 'bold'

  data$adj_pval <- -log(data$adj_pval, 10)
	suby <- data[!duplicated(data$term), ]
	if (is.numeric(nsub) == T){		
		suby <- suby[1:nsub, ]
	}else{
	  if (strsplit(nsub[1], ':')[[1]][1] == 'GO'){
	    suby <- suby[suby$ID%in%nsub, ]
	  }else{
	    suby <- suby[suby$term%in%nsub, ]
	  }
		nsub <- length(nsub)}
	N <- dim(suby)[1]
	r_pval <- round(range(suby$adj_pval), 0) + c(-2, 2)
	ymax <- c()
	for (i in 1:length(suby$adj_pval)){
		val <- (suby$adj_pval[i] - r_pval[1]) / (r_pval[2] - r_pval[1])
		ymax <- c(ymax, val)}
	df <- data.frame(x = seq(0, 10 - (10 / N), length = N), xmax = rep(10 / N - 0.2, N), y1 = rep(rad1, N), y2 = rep(rad2, N), ymax = ymax, zscore = suby$zscore, ID = suby$ID)
	scount <- data[!duplicated(data$term), which(colnames(data) == 'count')][1:nsub]
	idx_term <- which(!duplicated(data$term) == T)
	xm <- c(); logs <- c()
	for (sc in 1:length(scount)){
		idx <- c(idx_term[sc], idx_term[sc + 1] - 1)
		val <- runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06))
	 	xm <- c(xm, val)
		r_logFC <- round(range(data$logFC[idx[1]:idx[2]]), 0) + c(-1, 1)
			for (lfc in idx[1]:idx[2]){
				val <- (data$logFC[lfc] - r_logFC[1]) / (r_logFC[2] - r_logFC[1])
				logs <- c(logs, val)}
	}
	cols <- c()
	for (ys in 1:length(logs)) cols <- c(cols, ifelse(data$logFC[ys] > 0, 'upregulated', 'downregulated'))
	dfp <- data.frame(logx = xm, logy = logs, logFC = factor(cols), logy2 = rep(rad2, length(logs)))
	c <-	ggplot()+
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y1, ymax = y1 + ymax, fill = zscore), colour = 'black') +
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2, ymax = y2 + 1), fill = 'gray70') +
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.5, ymax = y2 + 0.5), colour = 'white') +
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.25, ymax = y2 + 0.25), colour = 'white') +
    geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.75, ymax = y2 + 0.75), colour = 'white') +
    geom_text(data = df, aes(x = x + (xmax / 2), y = y2 + 1.3, label = ID, angle = 360 - (x = x + (xmax / 2)) / (10 / 360)), size = label.size, fontface = label.fontface) +
    coord_polar() +
    labs(title = title) +
    ylim(1, rad2 + 1.6) +
    xlim(0, 10) +
    theme_blank +
    scale_fill_gradient2('z-score', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), breaks = c(min(df$zscore), max(df$zscore)),labels = c('decreasing', 'increasing')) +
    theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal') +	
    geom_point(data = dfp, aes(x = logx, y = logy2 + logy), pch = 21, fill = 'transparent', color = 'black', size = 3)+
    geom_point(data = dfp, aes(x = logx, y = logy2 + logy, color = logFC), size = 2.5)+
    scale_colour_manual(values = lfc.col, guide = guide_legend(title.position = "top", title.hjust = 0.5))		

	if (table.legend){
		table <- draw_table(suby, size = 4, col = 'black')
		c <- c + theme(plot.margin = unit(c(0.1, 0.1, 0.1, -1.5), 'cm'))
		par(mar = c(0.1, 0.1, 0.1, 0.1))
		grid.arrange(c, table, heights = c(2, 1), widths = c(1.8, 1), ncol = 2)
	}else{
		c + theme(plot.background = element_rect(fill = 'aliceblue'), panel.background = element_rect(fill = 'white'))
	}
}	

#' 
#' @name GOChord
#' @title Displays the relationship between genes and terms.
#' @description The GOChord function generates a circularly composited overview 
#'   of selected/specific genes and their assigned processes or terms. More 
#'   generally, it joins genes and processes via ribbons in an intersection-like
#'   graph. The input can be generated with the \code{\link{chord_dat}} 
#'   function.
#' @param data The matrix represents the binary relation (1= is related to, 0= 
#'   is not related to) between a set of genes (rows) and processes (columns); a
#'   column for the logFC of the genes is optional
#' @param space The space between the chord segments of the plot
#' @param gene.order A character vector defining the order of the displayed gene
#'   labels
#' @param gene.size The size of the gene labels
#' @param gene.space The space between the gene labels and the segement of the 
#'   logFC
#' @param nlfc Defines the number of logFC columns (default=1)
#' @param lfc.col The fill color for the logFC specified in the following form: 
#'   c(color for low values, color for the mid point, color for the high values)
#' @param lfc.min Specifies the minimium value of the logFC scale (default = -3)
#' @param lfc.max Specifies the maximum value of the logFC scale (default = 3)
#' @param ribbon.col The background color of the ribbons
#' @param border.size Defines the size of the ribbon borders
#' @param process.label The size of the legend entries
#' @param title The title (on top) of the plot
#' @details The \code{gene.order} argument has three possible options: "logFC", 
#'   "alphabetical", "none". Actually, the options are quite self- explanatory.
#'   
#'   Another argument which needs a bit more explanation is \code{nlfc}.
#'   Differential expression analysis can be performed for multiple conditions 
#'   and/or batches. Therefore, the data frame contains more than one logFC
#'   value per gene. To adjust to this situation the \code{nlfc} argument is
#'   used. It is a numeric value and it defines the number of logFC columns
#'   within the binary membership matrix. The default is "1" assuming that most
#'   of the time only one contrast is considered.
#' @seealso \code{\link{chord_dat}}
#' @import ggplot2
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Generating the binary matrix
#' chord<-chord_dat(circ,EC$genes,EC$process)
#' 
#' #Creating the chord plot
#' GOChord(chord)
#' 
#' #Creating the chord plot genes ordered by logFC and a different logFC color scale
#' GOChord(chord,space=0.02,gene.order='logFC',lfc.col=c('red','black','cyan'))
#' }
#' @export

GOChord <- function(data, title, space, gene.order, gene.size, gene.space, nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, process.label){
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  
  if (missing(title)) title <- ''
  if (missing(space)) space = 0
  if (missing(gene.order)) gene.order <- 'none'
  if (missing(gene.size)) gene.size <- 3
  if (missing(gene.space)) gene.space <- 0.2
  if (missing(lfc.col)) lfc.col <- c('brown1', 'azure', 'cornflowerblue')
  if (missing(lfc.min)) lfc.min <- -3
  if (missing(lfc.max)) lfc.max <- 3
  if (missing(ribbon.col)) colRib <- rainbow(Ncol - nlfc) else colRib <- ribbon.col
  if (missing(border.size)) border.size <- 0.5
  if (missing (process.label)) process.label <- 11
  
  if (gene.order == 'logFC') data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == 'alphabetical') data <- data[order(rownames(data)), ]
	if (sum(!is.na(match(colnames(data), 'logFC'))) > 0){
    if (nlfc == 1){
      cdata <- data[, 1:(Ncol - 1)]
		  lfc <-data[, Ncol]
    }else{
      cdata <- data[, 1:(Ncol - nlfc)]
      lfc <- data[, (Ncol - nlfc + 1)]
    }
	}else{
		cdata <- data
		lfc <- 0
	}
	nrib <- colSums(cdata)
	ngen <- rowSums(cdata)
	Ncol <- dim(cdata)[2]
	Nrow <- dim(cdata)[1]
	colRibb <- c()
	for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 202 * nrib[b]))
	r1 <- 1; r2 <- r1 + 0.1
	xmax <- c(); x <- 0
	for (r in 1:length(nrib)){
		perc <- nrib[r] / sum(nrib)
		xmax <- c(xmax, (pi * perc) - space)
		if (length(x) <= Ncol - 1) x <- c(x, x[r] + pi * perc)
	}
	xp <- c(); yp <- c()
	l <- 50
	for (s in 1:Ncol){
		xh <- seq(x[s], x[s] + xmax[s], length = l)
		xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), r2 * sin(x[s]))
		yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), r2 * cos(x[s]))
	}
	df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), each = 4 + 2 * l))
	xp <- c(); yp <- c(); logs <- NULL
	x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
	xmax2 <- rep(-pi / Nrow + space, length = Nrow)
	for (s in 1:Nrow){
		xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc == 1){
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), r2 * sin(x2[s]))
		  yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), r2 * cos(x2[s]))
    }else{
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc){
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
    }}}
  if(lfc[1] != 0){
    if (nlfc == 1){
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 2 * l))
    }else{
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc*Nrow)), each = 4 + 2 * l), logFC = rep(logs, each = 4 + 2 * l))  
      }
	}else{
	  df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l))
	}
	aseq <- seq(0, 180, length = length(x2)); angle <- c()
	for (o in aseq) if((o + 270) <= 360) angle <- c(angle, o + 270) else angle <- c(angle, o - 90)
	df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + xmax2/2),ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),labels = rownames(cdata), angle = angle)
	df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax / 2),ypro = (r1 + 0.15) * cos(x + xmax / 2), labels = colnames(cdata), stringsAsFactors = FALSE)
	cols <- rep(colRib, each = 4 + 2 * l)
	x.end <- c(); y.end <- c(); processID <- c()
	for (gs in 1:length(x2)){
		val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 1)
		pros <- which((cdata[gs, ] != 0) == T)
		for (v in 1:(length(val) - 1)){
			x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
			y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
			processID <- c(processID, rep(pros[v], 2))
		}
	}
	df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
	df_bezier <- df_bezier[order(df_bezier$processID,-df_bezier$y.end),]
	x.start <- c(); y.start <- c()
	for (rs in 1:length(x)){
		val<-seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 1)
		for (v in 1:(length(val) - 1)){
			x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
			y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
		}
	}	
	df_bezier$x.start <- x.start
	df_bezier$y.start <- y.start
	df_path <- bezier(df_bezier, colRib)
  tmp <- sapply(df_genes$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
  logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
  df_genes$logFC <- logFC
	
g<- ggplot() +
  geom_polygon(data = df_process, aes(x, y, group=id), fill='gray70', inherit.aes = F,color='black') +
  geom_polygon(data = df_process, aes(x, y, group=id), fill=cols, inherit.aes = F,alpha=0.6,color='black') +	
  geom_point(aes(x = xpro, y = ypro, size = factor(labels, levels = labels), shape = NA), data = df_texp) +
  guides(size = guide_legend("GO Terms", ncol = 4, byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
  theme(legend.text = element_text(size = process.label)) +
  geom_text(aes(xgen, ygen, label = labels, angle = angle), data = df_texg, size = gene.size) +
  geom_polygon(aes(x = lx, y = ly, group = ID), data = df_path, fill = colRibb, color = 'black', size = border.size, inherit.aes = F) +		
  labs(title = title) +
	theme_blank
    
  if (nlfc >= 1){
    g + geom_polygon(data = df_genes, aes(x, y, group = id, fill = logFC), inherit.aes = F, color = 'black') +
      scale_fill_gradient2('logFC', low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), 
                           breaks = c(min(df_genes$logFC), max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), round(max(df_genes$logFC)))) +
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }else{
    g + geom_polygon(data = df_genes, aes(x, y, group = id), fill = 'gray50', inherit.aes = F, color = 'black')+
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }
}

