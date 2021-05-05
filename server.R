require(dplyr)
require(Matrix)
load("generated_files/fracmat.rda")
load("generated_files/anno_expr_mat.rda")
allgroups_clusts=unique(annomat$cluster_label[order(annomat$cluster_id)])
set_plot_dimensions <- function(width_choice, height_choice) {
  options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), ylabval='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=ylabval)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

server <- function(input, output) {

  output$dotPlot <- renderPlot({
    input$go
    ###generate plotting matrix####
    inc=1/6
    genes=isolate(getgenelist())
    clusts=isolate(clustsetInput())
    if (!(tolower(clusts[1]) %in% c("all","glia","neuron"))) {
      clusts=annomat$cluster_label[match(clusts,annomat$cluster_id)]
    } else {
      if (clusts[1]=="all") {
        clusts=allgroups_clusts
      }
      if (tolower(clusts[1])=="neuron") {
        clusts=grep("GABA|Glut|Ex|Inh",allgroups_clusts,val=T)
      }
      if (tolower(clusts[1])=="glia") {
        clusts=grep("GABA|Glut|Ex|Inh",allgroups_clusts,val=T,invert=T)
      }
    }
    scaletype=isolate(scaletypeInput())
    genes=rownames(normdat1_tpm)[match(tolower(genes),tolower(rownames(normdat1_tpm)))]
    keepcols=which(annomat$cluster_label %in% clusts)
    groups=clusts
    testvec=annomat$cluster_label
    submat=matrix(NA,nrow=length(genes),ncol=length(groups)+1)
    colnames(submat)=c(groups,"50% expression size")
    rownames(submat)=genes
    fracmat=submat
    numcells=c()
    keepinds=list()
    for (jj in groups) {
      keepinds[[jj]]=annomat$sample_name[intersect(which(testvec==jj),keepcols)]
      numcells=c(numcells,length(keepinds[[jj]]))
    }
    names(numcells)=groups
    for (ii in 1:length(genes)) {
      allvals=normdat1_tpm[genes[ii],keepcols]
      names(allvals)=colnames(normdat1_tpm)[keepcols]
      if (scaletype=="zscore") {
        allvals2=c(scale(allvals))
      }
      if (scaletype=="logcpm") {
        allvals2=c(log10(allvals+1))
      }
      if (scaletype=="zscorelogcpm") {
        allvals2=c(scale(log10(allvals+1)))
      }
      names(allvals2)=names(allvals)
      for (jj in groups) {
        submat[genes[ii],jj]=mean(allvals2[keepinds[[jj]]]) 
        fracmat[genes[ii],jj]=sum(allvals[keepinds[[jj]]]>0)/numcells[jj]
      }
    }
    fracmat[,ncol(fracmat)]=0.5
    xpos=rep(1:ncol(submat),each=nrow(submat))
    ypos=rep(1:nrow(submat),ncol(submat))
    colinds=c(submat)
    if (scaletype=="logcpm") {
      colinds[colinds>5]=5
      colinds=round(colinds*200/5)+1
      colvec=colorRampPalette(c("white","yellow","red"))
      colvec=colvec(201)
      colplot=colvec[colinds]
    } else {
      colinds[colinds>2]=2
      colinds[colinds< -2]=-2
      colinds=round((colinds+2)*200/(4))+1
      colvec=colorRampPalette(c("blue","white","red"))
      colvec=colvec(201)
      colplot=colvec[colinds]
    }
    sizemat=c(fracmat)
    dfx = data.frame(x=c(xpos), y=max(ypos)-c(ypos)+1, sizeval=sizemat,colsplot=colplot)
    dfx = dfx[which(!is.na(dfx$sizeval)),]
    dfx = dfx[dfx$sizeval>0,]
    set_plot_dimensions(width_choice = ncol(fracmat),height_choice = nrow(fracmat))
    layout(matrix(c(rep(1,14),3,2), 2, 8, byrow = FALSE))
    par("mar"=c(12,7,14,0))
    plot(c(1,ncol(submat)),c(1,nrow(submat)),pch='',xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
    #submat[is.na(submat)]=0
    for (ii in 1:nrow(submat)) {
      lines(c(range(which(!is.na(submat[ii,])))),rep(max(ypos)-ii,2),col="grey")
    }
    for (jj in 1:(ncol(submat)-1)) {
      upval=max(dfx$y[dfx$x==jj])
      lines(c(jj,jj),c(upval,0),col="grey")
    }
    with(dfx, symbols(x=x, y=y, circles=sqrt(sizeval), inches=inc, ann=F, bg=as.character(colsplot), fg="black", xlab=colnames(submat),add=T,xlim=c(1,ncol(submat)),ylim=c(1,nrow(submat)),xaxt='n',yaxt='n'))
    axis(3, at=1:ncol(submat),labels=c(colnames(submat)[-ncol(submat)],""),las=2,cex.axis=1.2)
    axis(1, at=1:ncol(submat),labels=c(paste0("n = ",numcells),colnames(submat)[ncol(submat)]),las=2,cex.axis=1.2)
    axis(2, at=1:nrow(submat),label=rev(genes),las=2,cex.axis=1)
    if (scaletype=="zscore") {
      color.bar(colvec,-2,2,ylabval="Mean z-score (over selected cells)")
      }
    if (scaletype=="logcpm") {
      color.bar(colvec,0,5,ylabval="Log10(CPM+1)")
    }
    if (scaletype=="zscorelogcpm") {
      color.bar(colvec,-2,2,ylabval = "Mean z-score of Log10(CPM+1) (over selected cells)")
    }
  })
  
  genesetInput <- reactive({
    strsplit(input$genelist,",| |;|, ")[[1]]
  })
  
  clustsetInput <- reactive({
    strsplit(input$clustlist,",| |;|, ")[[1]]
  })
  
  scaletypeInput <- reactive({
    input$scaletyp
  })

  modevalInput <- reactive({
    input$modeval
  })
  
  clustonInput <- reactive({
   strsplit(input$clustonlist,",| |;|, ")[[1]] 
  })
  
  clustoffInput <- reactive({
    strsplit(input$clustofflist,",| |;|, ")[[1]] 
  })
  
  getgenelist <- reactive({
    input$go
    modeval=isolate(modevalInput())
    if (modeval=="findgenesval") {
      onclust=isolate(clustonInput())
      if (!(tolower(onclust[1]) %in% c("all","glia","neuron"))) {
        onclust=annomat$cluster_label[match(onclust,annomat$cluster_id)]
      } else {
        if (tolower(onclust[1])=="all") {
          onclust=allgroups_clusts
        }
        if (tolower(onclust[1])=="neuron") {
          onclust=grep("GABA|Glut|Ex|In",allgroups_clusts,val=T)
        }
        if (tolower(onclust[1])=="glia") {
          onclust=grep("GABA|Glut|Ex|In",allgroups_clusts,val=T,invert=T)
        }
      }
    
      offclust=isolate(clustoffInput())
      if (!(tolower(offclust[1]) %in% c("all","glia","neuron"))) {
        offclust=annomat$cluster_label[match(offclust,annomat$cluster_id)]
      } else {
        if (tolower(offclust[1])=="all") {
          offclust=allgroups_clusts
        }
        if (tolower(offclust[1])=="neuron") {
          offclust=grep("GABA|Glut|Ex|In",allgroups_clusts,val=T)
        }
        if (tolower(offclust[1])=="glia") {
          offclust=grep("GABA|Glut|Ex|In",allgroups_clusts,val=T,invert=T)
        }
      }
      
      offclust=setdiff(offclust,onclust)
      if (length(offclust)==0) {
        offclust=setdiff(allgroups_clusts,onclust)
      }
      numgenes=isolate(numgeneInput())
      if (length(onclust)==1) {fracval1=fracmat[,onclust]} else {fracval1=apply(fracmat[,onclust],1,mean)}
      if (length(offclust)==1) {fracval2=fracmat[,offclust]} else {fracval2=apply(fracmat[,offclust],1,mean)}
      diffval=fracval1-fracval2
      diffval=diffval[order(-diffval)]
      diffval=diffval[1:max(1,min(numgenes,length(which(diffval>0))))]
      genes=names(diffval)
    } else {
      genes=isolate(genesetInput())
      genes=genes[tolower(genes) %in% tolower(rownames(normdat1_tpm))]
    }
    return(genes)
  })
  
  numgeneInput <- reactive({
    as.numeric(input$numgeneval)
  })
  
  
  output$caption <- renderText({
    printlist=setdiff(tolower(genesetInput()),tolower(rownames(normdat1_tpm)))
    if (length(printlist)>0)
    paste0(paste(printlist,collapse=",")," not found","\n")
  })
  
  output$caption2 <- renderText({
    input$go
    paste0(isolate(getgenelist()),collapse=",")
  })
}

