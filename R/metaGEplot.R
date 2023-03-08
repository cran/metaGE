#' @import utils
## quiets concerns of R CMD check re: the .'s that appear in pipelines
#if(getRversion() >= "2.15.1")  utils::globalVariables(c(".",">"))
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


#' Drawing a QQplot
#'
#' The function MakeQQplot displays the QQplot of the -log10(pvalues).
#' @param Pvalues A vector containing pvalues.
#' @param Name A name of the corresponding test. (optional)
#' @param Xrange A range for the x axis. (optional)
#' @param Yrange A range for the y axis. (optional)

MakeQQplot <- function(Pvalues,Name=NULL,Xrange=NULL,Yrange=NULL){
  n=sum(!is.na(Pvalues))
  mLog10Pval <- Pvalues %>%
    sort %>%
    log10 %>%
    `*`(-1)
  mLog10Pval.H0 <- -log10((1:n)/n)
  if(is.null(Xrange)){
    Xrange=range(mLog10Pval.H0)
  }
  if(is.null(Yrange)){
    Xrange=range(mLog10Pval)
  }
  plot(mLog10Pval.H0,mLog10Pval,
       xlab='Expected -log10(pval)',ylab='Observed -log10(pval)',main=Name,
       xlim=Xrange,ylim=Yrange,
       type='l',lwd=2.5,cex.axis=1.5,cex.main=2)
  abline(0,1,col=2,lwd=2.5)
  lines(-log10((1:n)/n),mLog10Pval,
        lwd=2.5)

}

#' Display visual checks of pvalues.
#'
#' The function metaGE.pvalplot displays the pvalue distribution and the QQplot of the -log10(pvalues).
#' @param Pvalues A vector containing pvalues.
#' @param Main The main to display.(optional)
#' @return No return value, the plot is displayed in the active graphics window.
#' @export
#' @examples
#' # Import the data
#' data("metaData")
#'
#' # Compute the inter-environment correlation matrix
#' matCorr <- metaGE.cor(metaData, Threshold = 0.8)
#'
#' # Fit the Fixed Effect model
#' FeDF <- metaGE.fit(metaData, matCorr, Method = "Fe")
#'
#' # Check the pvalues
#' metaGE.pvalplot(Pvalues = FeDF$PVALUE, Main= "Pvalue Fe")

metaGE.pvalplot <- function(Pvalues, Main=''){
  # #Save the user graphical parameters
  # oldpar <- par(no.readonly = TRUE)
  # on.exit(par(oldpar))

  if(!is.numeric(Pvalues)){stop(paste0(Pvalues,' must be a numeric vector.'))}
  # par(mfrow=c(1,2))
  hist(Pvalues,100,main= Main,xlab='Pvalues', cex.axis=1.5,cex.main=2)
  MakeQQplot(Pvalues,Main)
}


#' Draw the heatmap to see markers effects across environments.
#'
#' The function metaGE.heatmap displays the heatmap of the zscores, the estimated marker effects or the pvalues of each markers (in rows) in each environments (in columns).
#' @param Data A dataset containing the zscores, the effects or the pvalues of each marker (in rows) in each environment (in columns), as obtained from metaGE.fit.
#' @param Prefix The prefix of the score to display in the heatmap: 'Z.' for the zscores, 'EFFECT.' for the effects and 'PVAL.' for the pvalues.('Z.' by default)
#' @param EnvGroups A dataset containing the names of the environments (in the first column) and the groups to which the environments belong (in the second column). (optional)
#' @param QTLsVarName The name of the column indicating to which QTL the marker belongs.  (optional)
#' @param RowOrder A boolean specifying whether to reorder the markers or not. (TRUE by default)
#' @param ColOrder A boolean specifying whether to reorder the environments or not. (TRUE by default)
#' @param ShowDendogram A boolean specifying wether to show the clustering of the rows and/or the columns. (FALSE by default)
#' @param Colors A vector of three colors corresponding to the
#' @return The heatmap
#' @export
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @examples
#' require(dplyr)
#' # Import the data
#' data("metaData")
#'
#' # Compute the inter-environment correlation matrix
#' matCorr <- metaGE.cor(metaData, Threshold = 0.8)
#'
#' # Fit the Fixed Effect model
#' FeDF <- metaGE.fit(metaData, matCorr, Method = "Fe")
#'
#' # Control the FDR (here Benjamini-Hochberg)
#' Alpha <- 0.05
#' Signif <-  FeDF$PVALUE %>% p.adjust(method = "BH") %>% `<`(Alpha) %>% which
#'
#' # Draw the z-scores heatmap of the significant markers
#'heatmap <- metaGE.heatmap(Data = FeDF[Signif,],
#'                          Prefix = "Z.")

metaGE.heatmap <- function(Data, Prefix='Z.', EnvGroups=NULL, QTLsVarName=NULL,RowOrder=TRUE, ColOrder = TRUE, ShowDendogram=FALSE, Colors=c("red","black","green")){

  ## Checkings

  ### Checkings on Data
  if((!is_tibble(Data))&(!is.data.frame(Data))){
    stop('Data should be either a data.frame or a tibble')
  }

  ### Checkings on Prefix
  NbScoreCols <- names(Data) %>%
    str_detect(pattern = Prefix) %>%
    sum
  if(NbScoreCols==0){
    stop(paste0('Could not find columns with prefix',Prefix,' to make the plot'))
  }else{
    ColNames <- names(Data) %>%
      .[str_which(.,pattern = Prefix)]
  }

  ### Checkings on Groups
  if(!is.null(EnvGroups)){
    if(!is.data.frame(EnvGroups)){
      stop('EnvGroups should be a dataframe.')
    }else{
      ColNames.sel <- EnvGroups[,1] %>%
        paste0(Prefix,.) %>%
        intersect(.,ColNames)
      ColumnSplit <- EnvGroups[,2]
      if(is_empty(ColNames.sel)){
        stop('Environments in EnvGroups do not match with column names in Data.')
      }
    }
  } else {
    ColNames.sel <- ColNames
    ColumnSplit <- NULL
  }


  ### Checkings on QTLs
  if(!is.null(QTLsVarName)){
    if(!(QTLsVarName %in% colnames(Data))){
      stop(paste0(QTLsVarName,' is not a column in Data.'))
    }else {
      QTLs_names <- unique(Data %>% pull(QTLsVarName)) %>% as.character()
      RowSplit <- Data %>% pull(QTLsVarName)
    }
  }else{
    ### Ordering the markers
    Data <-Data %>% arrange(.data$CHR, .data$POS)
    QTLs_names <- character(0)
    RowSplit <- NULL
  }


  ## Select columns if needed
  Dataheat <- Data %>% select(one_of(ColNames.sel))

  ## Build the HmMatrix
  HmMatrix <- Dataheat %>% as.matrix

  ## At least 2 markers ?
  if(nrow(HmMatrix)<2){
    stop('At least 2 rows required to draw a heatmap')
  }else{
    if(min(HmMatrix)<0 & max(HmMatrix)>0){
      color_fun <- circlize::colorRamp2(c(min(HmMatrix), 0, max(HmMatrix)), Colors)
    }else{
      color_fun <- circlize::colorRamp2(c(-max(abs(HmMatrix)), 0, max(abs(HmMatrix))), Colors)
    }
    HM <- ComplexHeatmap::Heatmap(HmMatrix, name="Z-score", col=color_fun,
                                  row_title = QTLs_names, row_title_rot = 0,
                                  cluster_rows = RowOrder, cluster_columns = ColOrder, show_row_dend = ShowDendogram,show_column_dend = ShowDendogram,
                                  row_split = RowSplit ,column_split = ColumnSplit,
                                  row_gap = unit(3, "mm"), column_gap =  unit(2, "mm"),
                                  row_title_gp = grid::gpar(cex=2),column_names_gp = grid::gpar(cex=2),column_title_gp = grid::gpar(cex=1.5),
                                  heatmap_legend_param = list( grid_width = unit(1, "cm"),legend_height = unit(7, "cm"),labels_gp=grid::gpar(cex=1.5),title_gp=gpar(cex=2)))


    return(HM)
  }
}



#' Draw the Manhattan plot.
#'
#' The function metaGE.manhattan displays the Manhattan plot of the -log10(p-value) or the local score of each marker along the genome.
#' @param Data A dataset containing the columns: CHR, POS, MARKER and the variable to plot for each marker, as obtained from metaGE.fit.
#' @param VarName The name of the column containing the variable to plot, generally the p-value or a score.
#' @param Threshold A threshold in order to draw a "genome-wide significant" line.  (optional)
#' @param SigZones A dataset containing the significant zones to plot, as obtained from metaGE.lscore. Must have columns: CHR, POS, START, END. (optional)
#' @param Score A boolean. If FALSE, the -log10 of the variable is plotted, useful for plotting p-values. If TRUE, the raw values of the variable is plotted, useful for plotting scores. (FALSE by default)
#' @param AnnotateMarker A list of markers name to annotate in the plot. (optional)
#' @param Main The main to display. (optional)
#' @param col A character vector indicating which colors to alternate for different chromosomes.(c('grey', 'black') by default)
#' @param colSigZones A character indicating which color to plot the significant zones.('blue' by default)
#' @param Ylim Two numeric values, specifying the lower limit and the upper limit of the  y-axe scale. (optional)
#' @return The Manhattan plot
#' @export
#' @importFrom qqman manhattan
#' @importFrom ggrepel geom_label_repel
#' @examples
#' require(dplyr)
#' # Import the data
#' data("metaData")
#'
#' # Compute the inter-environment correlation matrix
#' matCorr <- metaGE.cor(metaData, Threshold = 0.8)
#'
#' # Fit the Fixed Effect model
#' FeDF <- metaGE.fit(metaData, matCorr, Method = "Fe")
#'
#' # Control the FDR (here Benjamini-Hochberg)
#' Alpha <- 0.05
#' Signif <-  FeDF$PVALUE %>% p.adjust(method = "BH") %>% `<`(Alpha) %>% which
#'
#' # Draw the corresponding manhattan plot
#' #PvalThresholdFe <- FeDF[Signif,]$PVALUE%>% max %>% max(.,0)
#' #manhattan_pval <- metaGE.manhattan(Data = FeDF,VarName = 'PVALUE',
#' #                             Threshold = PvalThresholdFe,
#' #                              Main = '-log10(Pval) alongside the chromosome Fe method')
#'
#'
#'# Compute the score local
#' xi <- 2
#' FeScore <- metaGE.lscore(FeDF,"PVALUE", xi)
#'
#' # Draw the corresponding manhattan plot
#' manhattan_lscore <- metaGE.manhattan(Data = FeScore$Data,VarName = 'SCORE',
#'                                      SigZones = FeScore$SigZones, Score = TRUE,
#'                                      Main = 'Local score alongside the chromosome Fe method')
#'
#'
#'
#'


metaGE.manhattan <- function(Data, VarName, Threshold=NULL,SigZones=NULL,Score = FALSE, AnnotateMarker=NULL, Main='', col=c("grey", "black"), colSigZones='blue', Ylim=NULL){

  ### Calculate cumulative position for all markers
  don <- Data %>%

    # Compute chromosome size
    group_by(.data$CHR) %>%
    summarise(chr_len=max(.data$POS)-min(.data$POS),
              chr_min = min(.data$POS)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(.data$chr_len))-.data$chr_len) %>%
    select(-.data$chr_len) %>%

    # Add this info to the initial dataset
    left_join(Data, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(.data$CHR, .data$POS) %>%
    mutate( POScum=.data$POS+.data$tot - .data$chr_min + 1)

  ### Annotate some markers
  don <- don %>% mutate(is_annotate="no")
  if(!is.null(AnnotateMarker)){
    don[which(don$MARKER %in% AnnotateMarker),]$is_annotate <- "yes"
  }


  nbCHR <- length(unique(Data$CHR))

  if(!is.null(SigZones)){
    ### Calculate cumulative position for significative zones
    sigzones <- Data %>%

      # Compute chromosome size
      group_by(.data$CHR) %>%
      summarise(chr_len=max(.data$POS)-min(.data$POS),
                chr_min = min(.data$POS)) %>%

      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(as.numeric(.data$chr_len))-.data$chr_len) %>%
      select(-.data$chr_len) %>%

      # Add this info to the initial dataset
      left_join(SigZones, ., by=c("Chr"="CHR")) %>%

      # Add a cumulative position of each SNP
      arrange(.data$Chr, .data$Start) %>%
      mutate( Startcum=.data$Start+.data$tot - .data$chr_min + 1,
              Endcum=.data$End+.data$tot- .data$chr_min +1)
  }


  ### Compute the legend of the x axis
  axisdf = don %>% group_by(.data$CHR) %>% summarize(center=(max(.data$POScum) + min(.data$POScum)) / 2 )

  ### Rename the variable to plot
  don <- don %>% select(.data$CHR,.data$POS, .data$is_annotate, .data$POScum, all_of(VarName)) %>% rename("VarToPlot"=all_of(VarName))

  if(!Score){
    ### Built the plot
    manhattan <- ggplot() +
      # Show all points
      geom_point(data=don, aes(x=.data$POScum, y=-log10(.data$VarToPlot),color=as.factor(.data$CHR)), alpha=0.8, size=2) +
      scale_color_manual(values = rep_len(col, nbCHR)) +

      # custom X axis:
      scale_x_continuous( labels = axisdf$CHR, breaks= axisdf$center ) + xlab("Chromosome") +
      # remove space between plot area and x axis
      ylab(paste0("-log10(",VarName,")")) +

      # Custom the theme:
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )

    # Add label using ggrepel to avoid overlapping
    if(!is.null(AnnotateMarker)){
      manhattan <- manhattan +  ggrepel::geom_label_repel( data=subset(don, don$is_annotate=="yes"), aes(x=.data$POScum, y= -log10(min(.data$VarToPlot,Ylim[2])),label=str_sub(paste0("QTL",.data$CHR,"_",round(.data$POS / 1e6 , digits = 2)),end=-2)),fontface = 'bold', size=12)
    }
    # Add a threshold
    if(!is.null(Threshold)){
      manhattan <- manhattan + geom_hline(yintercept=-log10(Threshold), color = "red",size=1.5)
    }
    # Add limit to the y axes
    if(is.null(Ylim)){
      Ylim <- c(min(-log10(don$VarToPlot)), max(-log10(don$VarToPlot)))
    }
  }else{
    ### Built the plot
    manhattan <- ggplot() +
      # Show all points
      geom_point(data=don, aes(x=.data$POScum, y=.data$VarToPlot,color=as.factor(.data$CHR)), alpha=0.8, size=2) +
      scale_color_manual(values = rep_len(col, nbCHR )) +

      # custom X axis:
      scale_x_continuous( labels = axisdf$CHR, breaks= axisdf$center ) + xlab("Chromosome") +
      # remove space between plot area and x axis
      ylab(VarName) +

      # Custom the theme:
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    # Add label using ggrepel to avoid overlapping
    if(!is.null(AnnotateMarker)){
      manhattan <- manhattan +  ggrepel::geom_label_repel( data=subset(don, don$is_annotate=="yes"), aes(x=.data$POScum, y=min(.data$VarToPlot,Ylim[2]),label=str_sub(paste0("QTL",.data$CHR,"_",round(.data$POS / 1e6 , digits = 2)),end=-2)),fontface = 'bold', size=12)
    }
    # Add a threshold
    if(!is.null(Threshold)){
      manhattan <- manhattan + geom_hline(yintercept=Threshold, color = "red",size=1.5)
    }
    # Add limit to the y axes
    if(is.null(Ylim)){
      Ylim <- c(min(don$VarToPlot), max(don$VarToPlot))
    }
    if(!is.null(SigZones)){
      # Compute chromosome size
      width_box <- Data %>% group_by(.data$CHR) %>%
        summarise(chr_len=max(.data$POS)-min(.data$POS)) %>%
        pull(.data$chr_len) %>% mean() * 0.05

      for (i in 1:nrow(sigzones)){
        manhattan <- manhattan + geom_rect(data=data.frame(xmin=max(sigzones$Endcum[i]-width_box,0),
                                                           xmax=min(sigzones$Endcum[i]+width_box,max(don$POScum)),
                                                           ymin=max(0,Ylim[1]),
                                                           ymax=min(sigzones$LocalScoreMax[i],Ylim[2])),
                                           aes(xmin=.data$xmin,xmax=.data$xmax,ymin=.data$ymin,ymax=.data$ymax),
                                           fill=colSigZones,alpha=0.2)
      }
    }
  }

  # Coef <- 2
  # TH <- theme(axis.text=element_text(size=10*Coef),
  #             axis.title=element_text(size=15*Coef),
  #             plot.title=element_text(size=15*Coef),
  #             strip.text.x = element_text(size = 20*Coef),
  #             legend.text=element_text(size=12*Coef),
  #             legend.title = element_text(size=15*Coef))

  return(manhattan+ ylim(Ylim))

}




#' Plot the z-score of a marker according to a covariate.
#'
#' The function metaGE.regplot displays the graph of the z-scores of a marker according to a covariate.
#' @param Data A dataset containing the columns: MARKER and the z-scores of each marker (in rows) in each environment (in columns), as obtained from metaGE.collect.
#' @param Covariate  A dataset containing the values of one or more covariates (in columns) in each environment (in rows).
#' @param EnvName The name of the column containing the names of the environment in the Covariate dataset.
#' @param MarkerName The name of the marker.
#' @param VarName The name of the column containing the covariable to plot.
#' @param Zscore A boolean. If FALSE, the estimated marker effects is plotted. If TRUE, the z-scores of the marker is plotted. (FALSE by default)
#' @param aesCol The name of the column containing a qualitative covariable to specify the color of the points. (optional)
#' @param Main The main to display.(optional)
#' @return The plot
#' @export
#' @import ggplot2 viridis
#' @importFrom ggrepel geom_text_repel
#' @examples
#' data("metaData")
#' data("envDesc")
#' metaGE.regplot(Data = metaData, Covariate = envDesc, EnvName = "ShortName",
#'                MarkerName = "AX-91369217", VarName = "Tnight.mean", aesCol = "Classification")
#'

metaGE.regplot <- function(Data, Covariate, EnvName, MarkerName, VarName, Zscore=FALSE, aesCol=NULL, Main=''){

  #Rename the variable
  Covariate <- Covariate %>%  rename("ENV" = all_of(EnvName))
  Covariate <- Covariate %>%  rename("COVARIABLE" = all_of(VarName))
  if(!is.null(aesCol)){ Covariate <- Covariate %>%  rename("AesCol" = all_of(aesCol))}

  #Selection of the marker
  SNP_df <- Data %>% filter(.data$MARKER == MarkerName)

  #### Jointure with the covariate
  if(!Zscore){
    SNP_df <- SNP_df %>% select(.data$MARKER, contains("EFFECT.")) %>% rename_with(~str_remove(.x,pattern = 'EFFECT.'))
  }
  if(Zscore){
    SNP_df <- SNP_df %>% select(.data$MARKER, contains("Z.")) %>% rename_with(~str_remove(.x,pattern = 'Z.'))
  }

  #Reshape in long
  SNP_df <- SNP_df %>% gather(key = "ENV", value = "EFFECT", -.data$MARKER)

  #Add to snp effect table
  SNP_df <- SNP_df %>% left_join(Covariate, ., by=c("ENV"="ENV"))

  #### Plot
  if(!is.null(aesCol)){
    plott <- ggplot(SNP_df, aes(x = .data$COVARIABLE, y = .data$EFFECT, color= factor(.data$AesCol))) +
      scale_color_viridis(discrete = TRUE, option = "D")+
      labs(color = aesCol)
  }else{
    plott <- ggplot(SNP_df, aes(x = .data$COVARIABLE, y = .data$EFFECT))
  }

  plott <- plott +
    geom_point(size = 5)  +
    geom_hline(yintercept = 0 , linetype = 2) +
    ggrepel::geom_text_repel(aes(label = .data$ENV), size=7,colour="black")+
    xlab(VarName) + ylab("Marker effect") +
    ggtitle(Main)

  if(Zscore){plott <- plott + ylab("Z-score") }

  return(plott)

}




