#BCyto----
#rm(list = ls())

#STEPS:
#upload the R file to the github and try to install from it (using devtools)
#write things in the Help and About buttons
#re-install everything in RStudio and try to run the app just with the libraries that I'm invoking
#test it in R without RStudio
#check dependencies



#NEXT: need to add all "package::" to the third party functions... (then I can remove all library calls but devtools) [go on the STOP tags to see where I stopped]
#BUG: now I need to write for the packages I'm using to be installed in case they are not!





#add:
#spillover(




#library(devtools)



#' Runs BCyto app
#' @export
runBCyto = function(){


library(shiny)
library(shinyjs) #shinyjs::hide()/show() through useShinyjs()
library(flowWorkspace) #reads FCS files and provides gatingsets
library(flowCore) #draws gate tools, such as flowCore::rectangleGate()
library(MASS) #kde2d()
library(Rgraphviz) #AgNode()
library(scales) #something for compgraphs()
library(rhandsontable) #comp matrix
library(autospill) #automated comp matrix calculation doi: 10.1038/s41467-021-23126-8
library(Rtsne)
library(dplyr)
library(shinyFiles) #folder selection





script = "$('#comp tr td').each(function(){
              var cellValue = $(this).text();
              if (cellValue == 100.0000000){
                $(this).css('background-color', 'lightgrey');
                $#(this).css('font-weight', 'normal') }
              if (cellValue == 0){
                $(this).css('font-weight', 'normal') }
              if (cellValue > 0 && cellValue < 10) {
                $(this).css('background-color', '#FFFFB2');
                $(this).css('font-weight', 'normal') }
              if (cellValue >= 10 && cellValue < 20) {
                $(this).css('background-color', '#FED976');
                $(this).css('font-weight', 'normal') }
              if (cellValue >= 20 && cellValue < 40) {
                $(this).css('background-color', '#FEB24C');
                $(this).css('font-weight', 'normal') }
              if (cellValue >= 40 && cellValue < 70) {
                $(this).css('background-color', '#FD8D3C');
                $(this).css('font-weight', 'normal') }
              if (cellValue >= 70 && cellValue < 100) {
                $(this).css('background-color', '#F03B20');
                $(this).css('font-weight', 'normal') }
              if (cellValue > 100) {
                $(this).css('background-color', '#BD0026');
                $(this).css('font-weight', 'normal') } })"
change_hook = "function(el,x) {
              var hot = this.hot;
              var cellChanges = [];
              var changefn = function(changes,source) {
              if (source === 'edit' || source === 'undo' || source === 'autofill' || source === 'paste') {
              row = changes[0][0];
              col = changes[0][1];
              oldval = changes[0][2];
              newval = changes[0][3];
              if (oldval !== newval) {
                var cell = hot.getCell(row, col);
                cell.style.color = 'red';
                cellChanges.push({'rowid':row, 'colid':col}); } } }
              var renderfn = function(isForced) {
              for(i = 0; i < cellChanges.length; i++){
              var rowIndex = cellChanges[i]['rowid'];
              var columnIndex = cellChanges[i]['colid'];
              var cell = hot.getCell(rowIndex, columnIndex);
              cell.style.color = 'red'; } }
              var loadfn = function(initialLoad) {
              for(i = 0; i < cellChanges.length; i++){
                    var rowIndex = cellChanges[i]['rowid'];
                    var columnIndex = cellChanges[i]['colid'];
                    var cell = hot.getCell(rowIndex, columnIndex);
                    cell.style.color = 'black'; }
              cellChanges = [] }
              hot.addHook('afterChange', changefn);
              hot.addHook('afterRender', renderfn);
              hot.addHook('afterLoadData', loadfn); }"
jscode = "shinyjs.disableTab = function(name) {
          var tab = $('.nav li a[data-value=' + name + ']');
          tab.bind('click.tab', function(e) {
          e.preventDefault();
          return false; });
          tab.addClass('disabled'); }
          shinyjs.enableTab = function(name) {
          var tab = $('.nav li a[data-value=' + name + ']');
          tab.unbind('click.tab');
          tab.removeClass('disabled'); }"
css = ".nav li a.disabled {
        color: #737373 !important;
        cursor: not-allowed !important; }"

fluochannelsfull <<- Y <<- ""
concat <<- gs <<- prolifchannel <<- prolifparent <<- NULL


#BUG: the app is allowing more than 4 gates in a thing...


initializing = function(filelist){
  setwd(filepath)
  cs = flowWorkspace::load_cytoset_from_fcs(filelist)
  uncompgs <<- flowWorkspace::GatingSet(cs)
  gs <<- flowWorkspace::gs_clone(uncompgs)
  compdfs <<- list(spillover(flowWorkspace::gs_get_cytoframe(gs, onlysampleIDs[1]))[[1]])
  names(compdfs) <<- "Cytometer-defined"
  rownames(compdfs[[1]]) <<- colnames(compdfs[[1]])
  flowCore::compensate(gs, compdfs[[1]])
  allchannels <<- colnames(cs)
  fluochannels <<- allchannels[grepl("SC", allchannels) != TRUE][allchannels[grepl("SC", allchannels) != TRUE] != "Time"]
  ch <<- as.list(allchannels)
  currentcustomaxis <<- list()
  for(i in seq(fluochannels)){
    currentcustomaxis[[i]] <<- c(4.42, 0, -100) }
  gs <<- transform(gs, transformerList(fluochannels, flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[1]][1], neg=currentcustomaxis[[1]][2], widthBasis=currentcustomaxis[[1]][3])))
  channelnames <<- flowWorkspace::cytoframe_to_flowFrame(cs[[length(cs)]])@parameters@data[,2]
  secondlist = list()
  for(i in seq(channelnames)){
    if(is.na(channelnames[i])){
      secondlist[[i]] = allchannels[i]
    } else{
      secondlist[[i]] = paste(allchannels[i], channelnames[i], sep = " :: ") } }
  names(ch) <<- unname(secondlist)
  X <<- allchannels[grep("FS", allchannels)[1]]
  Y <<- allchannels[grep("SS", allchannels)[1]]
  fluochannelsfull <<- names(ch[grepl("SC", ch) != TRUE][ch[grepl("SC", ch) != TRUE] != "Time"])
  results <<- data.frame(matrix(ncol = 1, nrow = length(filelist[onlysampleIDs])))
  rownames(results) <<- substr(filelist[onlysampleIDs], 1, nchar(filelist[onlysampleIDs]) - 4)
  colnames(results) <<- NA
  hTable <<- shiny::reactiveValues(data = as.data.frame(compdfs[[1]]*100))
  shiny::updateSelectInput(inputId = "sample", choices = c("", filelist[onlysampleIDs]), selected = filelist[onlysampleIDs][1])
  shiny::updateSelectInput(inputId = "previewsample", choices = c("", filelist), selected = filelist[onlysampleIDs][1])
  shiny::updateSelectInput(inputId = "bgsample", choices = c("", filelist[onlysampleIDs]))
  shiny::updateSelectInput(inputId = "prolifsamplebutt", choices = c("", filelist[onlysampleIDs]), selected = filelist[onlysampleIDs][1])
  shiny::updateSelectInput(inputId = "channelY", choices = c(ch), selected = Y)
  shiny::updateSelectInput(inputId = "channelX", choices = c(ch), selected = X)
  shiny::updateSelectInput(inputId = "ovchannelY", choices = c("", ch))
  shiny::updateSelectInput(inputId = "ovchannelX", choices = c("", ch))
  shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
  shinyjs::js$enableTab("plottab")
  shinyjs::js$enableTab("comptab") }
axisticks = function(axisID, type, customaxis){
  if(type == "log"){
    trans = flowWorkspace::flowjo_biexp(channelRange=4096, maxValue=262144, pos=customaxis[1], neg=customaxis[2], widthBasis=customaxis[3])
    majorticks = trans(c(-1000000, -100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000, 1000000))
    labels = expression(10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)
    minorticks = list()
    for(i in seq(labels)[-length(labels)]){
      if(i < 7){
        minorticks = append(minorticks, (majorticks[i+1] - majorticks[i])*(1 - log10(9:1))+majorticks[i])
      } else{
        minorticks = append(minorticks, (majorticks[i+1] - majorticks[i])*log10(1:9)+majorticks[i]) } }
    minorticks[length(minorticks) + 1] = (majorticks[7] - majorticks[6])*log10(1:9)[2] + minorticks[length(minorticks)][[1]] #it was 6 and 4
    axis(axisID, at = minorticks, lwd = 2, labels = FALSE, tck = -0.015)
    adaptedlabels = labels
    if(abs(majorticks[6] - majorticks[8]) < 600){ #hiding 10ˆ1
      adaptedlabels[6] = adaptedlabels[8] = "" }
    if(abs(majorticks[5] - majorticks[9]) < 600){ #hiding 10ˆ2
      adaptedlabels[5] = adaptedlabels[9] = "" }
    if(abs(majorticks[4] - majorticks[10]) < 600){ #hiding 10ˆ3
      adaptedlabels[4] = adaptedlabels[10] = "" }
    if(abs(majorticks[3] - majorticks[11]) < 600){ #hiding 10ˆ4
      adaptedlabels[3] = adaptedlabels[11] = "" }
    axis(axisID, at = majorticks, lwd = 2, labels = FALSE)
    axis(axisID, at = majorticks, lwd = 0, cex.axis = 1.1, line = -0.3, labels = adaptedlabels, las = 1) #labels
  } else if(type == "cont"){
    axis(axisID, at = seq(0, 250000, by = 50000), lwd = 2, labels = FALSE)
    axis(axisID, at = seq(0, 260000, by = 10000), lwd = 1, labels = FALSE, tck = -0.015)
    axis(axisID, at = seq(0, 250000, by = 50000), lwd = 0, cex.axis = 1.1, line = -0.3, labels = c("0", "50K", "100K", "150K", "200K", "250K"), las =1) #labels
  } else if(type == "histo"){
    majorticks = axis(2, labels = FALSE, lwd = 0)
    axis(2, at = majorticks, lwd = 2, labels = FALSE)
    axis(2, at = seq(0, max(majorticks)*1.2, by = majorticks[1]+majorticks[2]/5), lwd = 1, labels = FALSE, tck = -0.015) #minorticks
    shortlabels = NULL
    for(i in seq(majorticks)){
      if(majorticks[i] >= 1000){
        n = majorticks[i]
        n2 = n - as.numeric(substr(n, nchar(n)-2, nchar(n)))
        diminishednumber = n - n2
        if(diminishednumber/1000 == 0){
          shortlabels[[i]] = paste0(n2/1000, "K")
        } else{
          shortlabels[[i]] = paste0(n2/1000, ".", diminishednumber/100, "K") }
      } else{
        shortlabels[[i]] = majorticks[i] } }
    axis(2, at = majorticks, lwd = 0, cex.axis = 1.1, line = -0.3, las = 1, label = as.character(majorticks))
  } else if(type == "Overlaid histogram"){
    majorticks = axis(2, labels = FALSE, lwd = 0)
    axis(2, at = majorticks, lwd = 2, labels = FALSE)
    axis(2, at = seq(0, 1.2, by = majorticks[1]+majorticks[2]/4), lwd = 1, labels = FALSE, tck = -0.015)
    axis(2, at = majorticks, lwd = 0, cex.axis = 1.1, line = -0.3, las = 1, label = (majorticks*100)) } }
tsneaxisticks = function(Xlim, Ylim){
  axis(1, at = seq(Xlim[1]*0.9, Xlim[2]*0.9, length.out = 6), lwd = 2, labels = FALSE)
  axis(1, at = seq(Xlim[1]*0.9, Xlim[2]*0.9, length.out = 21), lwd = 1, labels = FALSE, tck = -0.015)
  axis(1, at = seq(Xlim[1]*0.9, Xlim[2]*0.9, length.out = 6), lwd = 0, cex.axis = 1.1, line = -0.3, labels = as.character(seq(0, 100, by = 20)), las = 1)
  axis(2, at = seq(Ylim[1]*0.9, Ylim[2]*0.9, length.out = 6), lwd = 2, labels = FALSE)
  axis(2, at = seq(Ylim[1]*0.9, Ylim[2]*0.9, length.out = 21), lwd = 1, labels = FALSE, tck = -0.015)
  axis(2, at = seq(Ylim[1]*0.9, Ylim[2]*0.9, length.out = 6), lwd = 0, cex.axis = 1.1, line = -0.3, labels = as.character(seq(0, 100, by = 20)), las = 1)
  title(xlab = "t-SNE 1", ylab = "t-SNE 2", cex.lab = 1.5, line = 2.5, font.lab = 2) }
newplot = function(sampleID, X, Y, parent, type, showaxis, axisfont, bgdf = NULL, ancestry = FALSE){
  drawplot = NULL
  if(length(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[sampleID], parent)[[1]])@exprs) == 0){
    drawplot = FALSE
  } else{
    drawplot = TRUE }
  if(drawplot == TRUE){
    if(grepl("SC", X)){
      Xlim = c(0, 264000)
    } else{
      Xlim = c(0, 4100) }
    dataf <<- as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[sampleID], parent)[[1]])@exprs)
    if(type == "Contour" || type == "Pseudocolor + Contour"){
      if(nrow(dataf) > 5000){
        set.seed(1)
        sampleddf = dataf[sample.int(nrow(dataf), 5000),]
      } else{
        sampleddf = dataf } }
    dfX = dataf[,X]
    dfX[dfX < 0] = 0
    if(type != "Histogram"){
      dfY = dataf[,Y]
      dfY[dfY < 0] = 0
      if(grepl("SC", Y)){
        Ylim = c(0, 264000)
      } else{
        Ylim = c(0, 4100) }
      if(type != "Contour"){
        if(type != "Backgating"){
          if(type == "Density"){
            col = c("grey0", "grey20", "grey40", "grey60", "grey80", "grey100")
          } else{
            col = c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142") }
          density = suppressWarnings(grDevices::densCols(dfX, dfY, colramp = colorRampPalette(col)))
          plot(dfX, dfY, xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, col = density, xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
        } else{
          plot(dfX, dfY, xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, col = "gray", xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
          points(bgdf[,X], bgdf[,Y], pch = 20, cex = 0.5, lwd = 0, col ="red") } }
      if(type == "Contour" || type == "Pseudocolor + Contour"){
        z <<- MASS::kde2d(sampleddf[,X], sampleddf[,Y], n = 100)
        if(type == "Contour"){
          plot(dfX, dfY, xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, col = "black", xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
          for(i in contourLines(z, nlevels = 12)){
            if(i$level != 5e-08){
              polygon(i$x, i$y, col = "white", lwd = 1.2) } }
        } else{
          contour(z, drawlabels = FALSE,  xaxt = "n", yaxt = "n", add = TRUE, ann = FALSE, lwd = 1.5, nlevels = 12) } }
      if(grepl("SC", Y)){
        axisticks(2, "cont", currentcustomaxis[[which(fluochannels == Y)]])
        if(showaxis == TRUE){
          title(ylab = Y, cex.lab = 1+axisfont/10, line = 3.3, font.lab = 2) }
      } else{
        axisticks(2, "log", currentcustomaxis[[which(fluochannels == Y)]])
        if(showaxis == TRUE){
          title(ylab = names(ch)[which(colnames(gs) == Y)], cex.lab = 1+axisfont/10, line = 2.6, font.lab = 2) } }
    } else{
      if(length(dfX) > 1){
        referencehist <<- hist(dfX, breaks = 20, plot = FALSE)
        refcounts = referencehist$counts
        refdensity = referencehist$density
        multiplier = refcounts/refdensity
        densline <<- density(dfX)
        densline$y <<- densline$y*multiplier[1]
        maxdens <<- max(densline[2]$y)
        if(Y == "Prolif"){
          multiplier = 1.2
        } else{
          multiplier = 1.05 }
        plot(densline, xaxt = "n", yaxt = "n", ann = FALSE, xlim = Xlim, ylim = c(0, maxdens*multiplier), xaxs = "i", yaxs = "i")
        polygon(densline, col = "grey", border = "black", lwd = 1)
        axisticks(2, "histo", currentcustomaxis[[which(fluochannels == Y)]])
      } else{
        plot(0, xaxt = "n", yaxt = "n", ann = FALSE, xaxs = "i", yaxs = "i", cex = 0) }
      if(showaxis == TRUE){
        line = length(strsplit(as.character(max(refcounts)), "")[[1]])*0.5 + 1.5
        title(ylab = "Count", cex.lab = 1+axisfont/10, line = line, font.lab = 2) } }
    if(grepl("SC", X)){
      axisticks(1, "cont", currentcustomaxis[[which(fluochannels == X)]])
    } else{
      axisticks(1, "log", currentcustomaxis[[which(fluochannels == X)]]) }
    if(showaxis == TRUE){
      if(ancestry == TRUE){
        line = 2
      } else{
        line = 3 }
      title(xlab = names(ch)[which(colnames(gs) == X)], cex.lab = 1+axisfont/10, line = line, font.lab = 2) }
  } else{
    shinyjs::alert(paste0("The sample ", "'", filelist[sampleID[i]], "'", " has 0 events in the selected parent. Please choose another sample or change the parent.")) } }
overlay = function(sampleIDs, X, Y, parent, ovtype, ovtone, showaxis, axisfont){
  drawplot = NULL
  for(i in seq(sampleIDs)){
    if(length(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[sampleIDs[i]], parent)[[1]])@exprs) == 0){
      drawplot = FALSE
    } else{
      drawplot = TRUE } }
  if(drawplot == TRUE){
    par(mar = c(4, 6, 1, 25), lwd = 2) #bottom, left, top, and right
    if(grepl("SC", X)){
      Xlim = c(0, 264000)
    } else{
      Xlim = c(0, 4100) }
    dfs = list()
    samplenames <<- list()
    for(i in seq(sampleIDs)){
      dfs = append(dfs, list(as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[sampleIDs[i]], parent)[[1]])@exprs)))
      samplenames <<- append(samplenames, substr(filelist[sampleIDs[i]], 1, nchar(filelist[sampleIDs[i]]) - 4)) }
    #setting colors
    if(ovtype == "Overlaid histogram"){
      transparency = 0.6
      colorlist = c("grey" = "transparent", "transparent" = "black", "transparent" = "black")[1:length(sampleIDs)] #col/border
      lyt = c(0, 1, 2)
      lwd = c(1, 3, 2)
    } else if(ovtype == "Offset histogram"){
      transparency = 0.5
      histlist = list()
      Ylim = c(0, 0)
      offsets = c(0)
      if(ovtone == "Colorful"){
        precolors = c("grey", "red2", "green4", "blue3", "orange","black", "wheat3", "brown", "orchid3", "salmon")
      } else{
        preprecolors = c("grey0", "grey10", "grey20", "grey30", "grey40", "grey50", "grey60", "grey70", "grey80", "grey90")
        precolors = c()
        precolormultiplier = length(preprecolors)/length(sampleIDs)
        for(i in seq(sampleIDs)){
          if(i == 1){
            precolors[i] = preprecolors[i]
          } else if(i != length(sampleIDs)){
            precolors[i] = preprecolors[precolormultiplier*i]
          } else{
            precolors[i] = preprecolors[length(preprecolors)] } } }
      alphacolors = c()
      for(i in seq(sampleIDs)){
        alphacolors[i] = rgb(col2rgb(precolors[i])[1], col2rgb(precolors[i])[2], col2rgb(precolors[i])[3], 255*transparency, maxColorValue = 255) }
      colorlist = rep("black", length(sampleIDs))
      names(colorlist) = alphacolors
      lyt = lwd = rep(1, length(sampleIDs))
    } else{
      if(ovtone == "Colorful"){
        colorlist = c("black", "red2")
      } else{
        colorlist = c("black", "grey") } }
    #setting data frames
    dfsX = list()
    for(i in seq(sampleIDs)){
      dfX = dfs[[i]][,X]
      dfX[dfX < 0] = 0
      dfsX = append(dfsX, list(dfX)) }
    if(ovtype == "Dot plot"){
      dfsY = list()
      for(i in seq(sampleIDs)){
        dfY = dfs[[i]][,Y]
        dfY[dfY < 0] = 0
        dfsY = append(dfsY, list(dfY)) }
      if(grepl("SC", Y)){
        Ylim = c(0, 264000)
      } else{
        Ylim = c(0, 4100) } }
    #making plots
    if(ovtype == "Dot plot"){
      dffinal = data.frame()
      lastcolum = c()
      for(i in seq(sampleIDs)){
        if(length(dffinal) == 0){
          dffinal = data.frame("X" = dfsX[[i]], "Y" = dfsY[[i]], "color" = colorlist[i])
        } else{
          dffinal = rbind(dffinal, data.frame("X" = dfsX[[i]], "Y" = dfsY[[i]], "color" = colorlist[i])) } }
      set.seed(2)
      shuffleddf = dffinal[sample(nrow(dffinal)), ]
      plot(shuffleddf[[1]], shuffleddf[[2]], xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, col = shuffleddf[[3]], xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")  }
    for(i in seq(sampleIDs)){
      if(ovtype == "Overlaid histogram" || ovtype == "Offset histogram"){
        densline = density(dfsX[[i]])
        densline$y = densline$y/max(densline$y)
        maxdens = max(densline$y)
        if(ovtype == "Overlaid histogram"){
          Ylim = c(0, maxdens*1.05)
          if(i == 1){
            plot(densline, col = "transparent", xaxt = "n", yaxt = "n", ann = FALSE, xlim = Xlim, ylim = c(0, maxdens*1.05), xaxs = "i", yaxs = "i") }
          polygon(densline, col = names(colorlist[i]), border = colorlist[i], lwd = lwd[i], lty = lyt[i])
        } else{
          if(i != 1){
            offsets[i] = Ylim[2] }
          Ylim[2] = Ylim[2] + maxdens*0.7 } } }
    if(ovtype == "Offset histogram"){
      Ylim[2] = Ylim[2] + maxdens*0.7/2*1.1
      plot(0, xaxt = "n", yaxt = "n", ann = FALSE, xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i", cex = 0)
      for(i in rev(seq(sampleIDs))){
        densline = density(dfsX[[i]])
        densline$y = densline$y/max(densline$y) + offsets[i]
        polygon(densline, col = names(colorlist[i]), border = NA, lwd = 1.5, lty = 1)
        lines(densline, lwd = 1.5, lty = 1) } }
    #adding axis
    if(grepl("SC", X)){
      axisticks(1, "cont", NULL)
    } else{
      axisticks(1, "log", currentcustomaxis[[which(fluochannels == X)]]) }
    if(ovtype == "Dot plot"){
      if(grepl("SC", Y)){
        distance = 3.3
        axisticks(2, "cont", NULL)
      } else{
        distance = 2.8
        axisticks(2, "log", currentcustomaxis[[which(fluochannels == Y)]]) }
    } else if(ovtype == "Overlaid histogram"){
      axisticks(2, ovtype, NULL) }
    if(showaxis == TRUE){
      if(ovtype == "Overlaid histogram" || ovtype == "Offset histogram"){
        title(xlab = names(ch)[which(colnames(gs) == X)], cex.lab = 1+axisfont/10, line = 3, font.lab = 2)
        if(ovtype == "Overlaid histogram"){
          title(ylab = "% of max", cex.lab = 1+axisfont/10, line = 2.5, font.lab = 2) }
      } else{
        title(ylab = names(ch)[which(colnames(gs) == Y)], cex.lab = 1+axisfont/10, line = distance, font.lab = 2)
        title(xlab = names(ch)[which(colnames(gs) == X)], cex.lab = 1+axisfont/10, line = 3, font.lab = 2) } }
    if(ovtype == "Overlaid histogram"){
      pch = c(22, NA, NA)[1:length(sampleIDs)]
      lty = c(0, 1, 2)[1:length(sampleIDs)]
    } else{
      pch = 22
      lty = 0 }
    if(ovtype == "Overlaid histogram" || ovtype == "Offset histogram"){
      pt.bg = names(rev(colorlist))
    } else{
      pt.bg = rev(colorlist) }
    templegend = legend("topright", bty = "n", inset = c(-0.25, -0.025), pch = rev(pch), lty = rev(lty), lwd = 2, cex = 1.3,pt.cex = 4, y.intersp = 1.6, text.width = strwidth("100"), pt.bg = pt.bg, col = rev(colorlist), legend = rep(NA, length(samplenames)), xpd = TRUE)
    multiplier = 1.2
    text(strwidth(rev(samplenames))*multiplier + templegend$text$x*1.01, templegend$text$y, rev(samplenames), pos = 2, xpd = NA, font = 2, cex = multiplier)
  } else{
    shinyjs::alert(paste0("The sample ", "'", filelist[sampleIDs[i]], "'", " has 0 events in the selected parent. Please choose another sample or change the parent.")) } }
detectgate = function(sampleID, X, Y, parent, type, observer, showgatename, gatefont, ancestry = FALSE){
  poppaths <<- flowWorkspace::gs_get_pop_paths(gs)

  detectedgate <<- detectedgatefullpath <<-NULL
  detectedgateinfo <<- list()
  gatetypes = NULL
  if(length(poppaths) > 1){
    for(i in seq(poppaths)[-1]){
      gatesandpopsinfo <<- flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[sampleID][[1]]
      gatetypes[i] = class(gatesandpopsinfo)[1]
      chX = names(gatesandpopsinfo@parameters[1])
      if(type != "Histogram"){
        if(length(gatesandpopsinfo@parameters) > 1){
          chY = names(gatesandpopsinfo@parameters[2])
        } else{
          chY = "none" }
      } else{
        chY = "Histogram" }
      if(X == chX && Y == chY || X == chY && Y == chX || type == "Histogram" && length(gatesandpopsinfo@parameters) == 1 && X == chX){
        detectedgateID = which(poppaths == poppaths[i])
        detectedparent <<- flowWorkspace::gs_pop_get_count_fast(gs[sampleID], "freq")[detectedgateID - 1][,3][[1]]
        if(parent == detectedparent){
          detectedgate[i-1] = gatesandpopsinfo@filterId
          detectedgatefullpath[i-1] <<- poppaths[i]
          detectedgateinfo[[i-1]] = flowWorkspace::gs_pop_get_gate(gs, detectedgate[i-1])[sampleID][[1]] } } }
    if(observer == "plot" || observer == ">= 4"){
      detectedgate <<- detectedgate[!is.na(detectedgate)]
      gatetypes = gatetypes[!is.na(gatetypes)]
      if(observer == "plot" && !is.null(detectedgate)){
        plotgate(sampleID, X, Y, type, detectedgate, gatetypes, showgatename, gatefont, ancestry) } } } }
plotgate = function(sampleID, X, Y, type, detectedgate, gatetypes, showgatename, gatefont, ancestry = FALSE){
  for(i in seq(detectedgate[!is.na(detectedgate)])){
    currentcoords = flowWorkspace::gs_pop_get_gate(gs, detectedgate[!is.na(detectedgate)][i])[[sampleID]]
    if(names(currentcoords@parameters[1]) == X){
      x = 1; y = 2
    } else{
      x = 2; y = 1 }
    axisXlim = par("usr")[2]
    axisYlim = par("usr")[4]
    if(type != "Histogram"){
      if(class(currentcoords)[1] == "polygonGate"){
        xleft = min(currentcoords@boundaries[,1])
        xright = max(currentcoords@boundaries[,1])
        ybottom = min(currentcoords@boundaries[,2])
        ytop = max(currentcoords@boundaries[,2])
        for(j in seq(nrow(currentcoords@boundaries))){
          if(j > 1){
            currentcoords@boundaries
            segments(currentcoords@boundaries[,1][j-1], currentcoords@boundaries[,2][j-1], currentcoords@boundaries[,1][j], currentcoords@boundaries[,2][j], lwd = 3) } }
      } else{
        if(currentcoords@min[[x]] == "-Inf"){
          xleft = -1000
        } else{
          xleft = currentcoords@min[[x]] }
        if(currentcoords@max[[x]] == "Inf"){
          xright = axisXlim*1.1
        } else{
          xright = currentcoords@max[[x]] }
        if(currentcoords@min[[y]] == "-Inf"){
          ybottom = -1000
        } else{
          ybottom = currentcoords@min[[y]] }
        if(currentcoords@max[[y]] == "Inf"){
          ytop = axisYlim*1.1
        } else{
          ytop = currentcoords@max[[y]] } }
      if((ybottom+ytop)/2 < max(axis(y, labels = FALSE, lwd = 0))/2){
        if(showgatename == TRUE){
          yjust = -0.5
        } else{
          yjust = -1.2 }
      } else{
        if(showgatename == TRUE){
          yjust = 1.5
        } else{
          yjust = 2.2 } }
      if(class(currentcoords)[1] == "rectangleGate"){
        rect(xleft = xleft, xright = xright, ybottom = ybottom,  ytop = ytop, lwd = 3) }
    } else{
      segments(currentcoords@min[[x]], axisYlim*0.75, currentcoords@max[[x]], axisYlim*0.75, lwd = 3)
      segments(currentcoords@min[[x]], axisYlim*0.72, currentcoords@min[[x]], axisYlim*0.78, lwd = 3)
      segments(currentcoords@max[[x]], axisYlim*0.72, currentcoords@max[[x]], axisYlim*0.78, lwd = 3) }
    gatedpopinfo = flowWorkspace::gs_pop_get_count_fast(gs[sampleID], subpopulations = detectedgate[!is.na(detectedgate)][i])
    frequency = paste0(sprintf("%.1f", gatedpopinfo[,4][[1]]*100/gatedpopinfo[,5][[1]]), "%")
    if(length(detectedgate[!is.na(detectedgate)]) == 4){
      if(ancestry == TRUE){
        QXcoords = c(0.2, 0.8, 0.8, 0.2)
        QYcoords = c(1.1, 1.1, 0.5, 0.5)
      } else{
        if(names(currentcoords@min)[1] == X){
          QXcoords = c(0.1, 0.9, 0.9, 0.1)
          QYcoords = c(1.05, 1.05, 0.25, 0.25)
        } else{
          QXcoords = c(0.9, 0.9, 0.1, 0.1)
          QYcoords = c(0.25, 1.05, 1.05, 0.25) } }
      if(showgatename == TRUE){
        if(detectedgate[!is.na(detectedgate)][i] == paste0("Q1: ", X, "- ", Y, "+") || detectedgate[!is.na(detectedgate)][i] == paste0("Q2: ", X, "+ ", Y, "+") || detectedgate[!is.na(detectedgate)][i] == paste0("Q3: ", X, "+ ", Y, "-") || detectedgate[!is.na(detectedgate)][i] == paste0("Q4: ", X, "- ", Y, "-")){
          title = substring(detectedgate[!is.na(detectedgate)][i], 1, 2)
        } else{
          title = detectedgate[!is.na(detectedgate)][i] }
      } else{
        title = NULL }
      legend(axisXlim*QXcoords[i], axisYlim*QYcoords[i], title = title, legend = frequency, cex = 1+gatefont/10, bg = rgb(1,1,1, 0.6), box.lwd = 0, x.intersp = -0.5, y.intersp = 0.8, text.font = 2, xjust = 0.5, yjust = 1.5)
    } else{
      if(showgatename == TRUE){
        title = detectedgate[!is.na(detectedgate)][i]
      } else{
        title = NULL }
      if(type != "Histogram"){
        legend((xleft+xright)/2, (ybottom+ytop)/2, title = title, legend = frequency, cex = 1+gatefont/10, bg = rgb(1,1,1, 0.6), box.lwd = 0, x.intersp = -0.5, y.intersp = 0.8, text.font = 2, xjust = 0.5, yjust = yjust)
      } else{
        legend((currentcoords@min[[x]]+currentcoords@max[[x]])/2, axisYlim*0.95, title = title, legend = frequency, cex = 1+gatefont/10, bg = rgb(1,1,1, 0.6), box.lwd = 0, x.intersp = -0.5, y.intersp = 0.8, text.font = 2, xjust = 0.5) } } } }
compgraphs = function(sampleID, newmatrix, overlay, inputpreview = NULL, parent = "root"){
  shiny::withProgress(message = "Please wait...", detail = "", value = 0, max = 100, {
    setwd(filepath)
    graphnumber = length(fluochannels) - 1
    par(oma = c(1, 2, 2, 1), mar = c(0, 0, 0, 0), mfrow = c(graphnumber, graphnumber), lwd = 2)
    cf = flowWorkspace::load_cytoframe_from_fcs(filelist[sampleID])
    originalmatrix = compdfs[[1]]
    rownames(originalmatrix) = colnames(originalmatrix)
    if(is.null(newmatrix)){
      flowCore::compensate(cf, originalmatrix)
    } else{
      flowCore::compensate(cf, newmatrix) }
    for(i in seq(fluochannels)){
      cf = transform(cf, transformerList(fluochannels[i], flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3]))) }
    dataf = as.data.frame(flowWorkspace::cytoframe_to_flowFrame(cf)@exprs)
    uncompdf = as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(uncompgsfortable[sampleID], parent)[[1]])@exprs)
    col = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
    col = alpha(col, alpha = 0.4)
    firstway = c()
    for(i in seq(graphnumber)){
      firstway[[i]] = as.numeric(newmatrix[i,])[-1:-i] }
    firstway = unlist(firstway)*100
    secondway = c()
    for(i in seq(graphnumber)){
      secondway[[i]] = as.numeric(newmatrix[,i])[-1:-i] }
    secondway = unlist(secondway)*100
    numberstoretrieve = list()
    for(i in seq(unlist(firstway))){
      numberstoretrieve[[i]] = c(unlist(firstway)[i], unlist(secondway)[i]) }
    ranges = c()
    for(i in 1:11){
      ranges[[i]] = c(9*i-9, 9*i) }
    lim = c(0, 4100)
    index = 0
    for(i in seq(graphnumber)){
      for(j in seq(graphnumber + 1 - i)){
        X = fluochannels[j+i]
        Y = fluochannels[i]
        set.seed(3)
        sampleddf = dataf[sample.int(nrow(dataf), 1000),]
        uncompsampleddf = uncompdf[sample.int(nrow(dataf), 1000),]
        index = index + 1
        plot(sampleddf[,X], sampleddf[,Y], xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, xlim = lim, ylim = lim, xaxs = "i", yaxs = "i")
        for(l in 1:11){
          if(max(numberstoretrieve[[index]]) == 0){
            bgcol = "white"
          } else{
            if(max(abs(numberstoretrieve[[index]])) > ranges[[l]][1] && max(abs(numberstoretrieve[[index]])) <= ranges[[l]][2] && max(abs(numberstoretrieve[[index]])) <= 81){
              bgcol = col[l]
            } else if((max(abs(numberstoretrieve[[index]])) > 81)){
              bgcol = col[9] } } }
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bgcol)
        points(sampleddf[,X], sampleddf[,Y], xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, xlim = lim, ylim = lim, xaxs = "i", yaxs = "i")
        if(overlay == TRUE){
          points(uncompsampleddf[,X], uncompsampleddf[,Y], col = "grey", xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = 0.5, lwd = 0, xlim = lim, ylim = lim, xaxs = "i", yaxs = "i") }
        if(length(fluochannels) >= 10){
          textsize = 0.6
        } else if(length(fluochannels) < 10 && length(fluochannels) > 4){
          textsize = rev(5:9*0.12)[length(fluochannels) - 4]
        } else if(length(fluochannels) <= 4){
          textsize = 1.2 }
        if(i == 1){
          mtext(X, side = 3, padj = -0.5, cex = textsize) }
        if(j == 1 || i > j && i == j - 1){
          mtext(Y, side = 2, padj = -0.5, cex = textsize) } }
      if(i != graphnumber){
        for(k in seq(i)){
          plot.new() } } } }) }




#UI----
ui <<- fluidPage(
  shinyjs::useShinyjs(),
  shinyjs::extendShinyjs(text = jscode, functions = c('disableTab','enableTab')),
  shinyjs::inlineCSS(css),
  tags$style("#about {border-color:white; font-size:12px}"),
  headerPanel(fluidRow(column(width = 3, div(icon("virus"), "BCyto") ),
                       column(width = 2, offset = 7, align = "right", actionButton(inputId = "about", label = "About BCyto"))), ),
  tags$style("#tabs {border-color:black; top:25px}"), tags$style(shiny::HTML(".tabbable > .nav > li > a {color:black}")), tags$style(shiny::HTML(".tabbable > .nav > li[class=active] > a {color:black; border-top-color:black; border-left-color:black; border-right-color:black}")), tags$style("#exportimagePlot {border-color:black}"), tags$style("#help {border-color:white}"),
  tags$style("#nextsample .fa{font-size: 11px}"), tags$style("#previoussample .fa{font-size: 11px}"), tags$style("#nextprolifsample .fa{font-size: 11px}"), tags$style("#prevprolifsample .fa{font-size: 11px}"),
  tabsetPanel(id = "tabs",
              ##File----
              tabPanel(div(icon("cloud-upload-alt"), "File"), value = "filetab",
                       br(),
                       fluidRow(column(width = 4, offset = 4,
                                       sidebarPanel(id = "Fileleft", align = "center",
                                                    shinyFiles::shinyDirButton("directory", "Load FCS files", icon = icon("file-upload"), "Please select a folder with your FCS files"),
                                                    br(), br(),
                                                    strong("or"),
                                                    br(), br(),
                                                    shinyFiles::shinyFilesButton("BCytofileload", "Load BCyto file", icon = icon("file-upload"), "Please select a BCyto file (it needs to be located in the same folder of the FCS files)", multiple = FALSE),
                                                    br(), br(),
                                                    strong("or"),
                                                    br(), br(),
                                                    actionButton(inputId = "loadtestdata", label = div(icon("file-download"), "Download test data")),
                                                    br(), br(), br(), br(), br(), br(),
                                                    actionButton(inputId = "savefile", label = div(icon("save"), "Save")),
                                                    br(), br(),
                                                    tags$style("#filehelp {border-color:white}"),
                                                    fluidRow(column(width = 3, offset = 9, actionButton(inputId = "filehelp", icon("question-circle"))) ),
                                                    width = 12), )),
              ),
              ##Plot----
              tabPanel(div(icon("chart-area"), "Plot"), value = "plottab",
                       br(),
                       tags$style(".well {background-color:white; border-color:black}"),
                       sidebarPanel(id = "Plotleft",
                                    selectInput(inputId = "sample", label = div(icon("vial"), "Sample"), choices = c("")),
                                    tags$style("#customizeaxisY {border-color:black}"),
                                    fluidRow(column(width = 8, selectInput(inputId = "channelY", label = div(icon("ruler-combined"), "Y axis"), choices = c(""))  ),
                                             column(width = 4, style = "margin-top: 25px;",  actionButton(inputId = "customizeaxisY", label = "Scale")) ),
                                    tags$style("#customizeaxisX {border-color:black}"),
                                    fluidRow(column(width = 8, selectInput(inputId = "channelX", label = div(icon("ruler-combined"), "X axis"), choices = c(""))  ),
                                             column(width = 4, style = "margin-top: 25px;",  actionButton(inputId = "customizeaxisX", label = "Scale")) ),
                                    selectInput(inputId = "type", label = div(icon("chart-area"), "Plot type"), choices = c("Pseudocolor", "Contour", "Density", "Pseudocolor + Contour", "Histogram")),
                                    tags$style("#displayoptmain {border-color:black}"),
                                    fluidRow(column(width = 4, actionButton(inputId = "displayoptmain", label = div(icon("bars"), "Display options"))  ),
                                             column(width = 2, offset = 5, actionButton(inputId = "help", icon("question-circle"))) ),
                                    width = 3),
                       mainPanel(id = "Plotmiddle",
                                 fluidRow(column(width = 3,
                                                 div(style = "display:inline-block; position:fixed",
                                                     actionButton(inputId = "previoussample", width = "40px", icon("chevron-left"), style="color: white; background-color:black; border-color:black"),
                                                     actionButton(inputId = "nextsample", width = "40px", icon("chevron-right"), style="color: white; background-color:black; border-color:black") ) ),
                                          column(width = 5, offset = 4, align = "left", id = "gatetools",
                                                 div(style = "display:inline-block; position:fixed",
                                                     actionButton(inputId = "rectangle", width = "40px", tags$i(class = "far fa-square"), style="color: white; background-color:black; border-color:black"),
                                                     actionButton(inputId = "polygon", width = "40px", icon("draw-polygon"), style="color: white; background-color:black; border-color:black"),
                                                     actionButton(inputId = "quadrant", width = "40px", icon("plus"), style="color: white; background-color:black; border-color:black"),
                                                     actionButton(inputId = "interval", width = "40px", icon("minus"), style="color: white; background-color:black; border-color:black"), ) ),
                                          column(width = 5, offset = 3,
                                                 div(style = "display:inline-block; position:fixed", textInput(inputId = "drawgatename", label = NULL, placeholder = "Please draw a gate", width = "180px"), ) ),
                                          column(width = 1, offset = 0,
                                                 div(style = "display:inline-block; position:fixed",
                                                     actionButton(inputId = "dotplotgateOK", label = "Ok", style="color: white; background-color:black; border-color:black"),
                                                     actionButton(inputId = "dotplotgateCANCEL", label = "Cancel", style="color: black; background-color:white; border-color:black") ) ), ),
                                 br(), br(),
                                 uiOutput("plotui"),
                                 tags$style("#saveMainPlot {border-color:black}"),
                                 div(style="display:inline-block", downloadButton(outputId = "saveMainPlot", label = "Export image")),
                                 div(style="display:inline-block;vertical-align:top; width: 20px;", shiny::HTML("<br>")),
                                 div(style="display:inline-block", textOutput("events")),
                                 width = 5),
                       sidebarPanel(id = "Plotright",
                                    strong(div(icon("sitemap"), "Parent")),
                                    plotOutput("hierarchy", click = "hierarchy_click", dblclick = "hierarchy_dblclick", height = 450),
                                    fluidRow(column(width = 6,
                                                    tags$style("#exportimageGates {border-color:black}"),
                                                    downloadButton(outputId = "exportimageGates", label = "Export image"), ),
                                             column(width = 5, offset = 1,
                                                    tags$style("#editgate {border-color:black}"),
                                                    actionButton(inputId = "editgate", label = "Edit", icon("fas fa-pen")), ), ),
                                    width = 4) ),
              ##Compensation----
              tabPanel(div(icon("table"), "Compensation"), value = "comptab",
                       br(),
                       sidebarPanel(
                         selectInput(inputId = "previewsample", label = div(icon("vial"), "Preview Sample"), choices = c("")),
                         checkboxInput(inputId = "showuncomp", label = strong("Overlay Uncompensated"), value = TRUE),
                         selectInput(inputId = "previewmatrix", label = div(icon("table"), "Preview Matrix"), choices = c("Cytometer-defined (applied)" = "Cytometer-defined")),
                         tags$style("#creatematrix {border-color:black}"),
                         actionButton(inputId = "creatematrix", label = "Create New Matrix", icon("plus-square")),
                         br(), br(),
                         tags$style("#savematrix {color: white; background-color:black; border-color:black}"),
                         actionButton(inputId = "savematrix", label = "Save New Matrix", icon("sd-card")),
                         br(), br(),
                         tags$style("#cancelmatrix {border-color:black}"),
                         actionButton(inputId = "cancelmatrix", label = "Cancel"),
                         br(), br(),
                         tags$style("#applymatrix {border-color:black}"),
                         actionButton(inputId = "applymatrix", label = "Apply to Samples", icon("share-square")),
                         width = 3),
                       mainPanel(
                         fluidRow(column(width = 12, rhandsontable::rHandsontableOutput("comp"))),
                         br(), br(),
                         fluidRow(column(width = 12, align = "center", plotOutput("compgraphs"))),
                         width = 9)
              ),
              ##Ancestry----
              tabPanel(div(icon("chart-area"), "Ancestry plots"), value = "ancestrytab",
                       br(),
                       sidebarPanel(
                         selectInput(inputId = "bgsample", label = div(icon("vial"), "Sample"), choices = c("")),
                         selectInput(inputId = "bgtype", label = div(icon("chart-area"), "Plot type"), choices = c("", "Pseudocolor", "Contour", "Density", "Backgating"), selected = ""),
                         selectInput(inputId = "bgpop", label = div(tags$i(class = "fas fa-clone"), shiny::HTML("Backgating <span style='color: red'>Population</span>")), choices = c("")),
                         width = 3),
                       mainPanel(
                         plotOutput("ancestry", height = 550),
                         tags$style("#exportimageAncestry {border-color:black}"),
                         downloadButton(outputId = "exportimageAncestry", label = "Export image"),
                         width = 9)
              ),
              ##Overlays----
              tabPanel(div(icon("chart-area"), "Overlays"), value = "overlaytab",
                       br(),
                       sidebarPanel(
                         selectInput(inputId = "ovtype", label = div(icon("chart-area"), "Plot type"), choices = c("", "Overlaid histogram", "Offset histogram", "Dot plot")),
                         selectInput(inputId = "ovtone", label = div(icon("chart-area"), "Tone"), choices = c("", "Colorful", "B&W")),
                         selectInput(inputId = "ovchannelY", label = div(icon("ruler-combined"), "Y axis"), choices = c("")),
                         selectInput(inputId = "ovchannelX", label = div(icon("ruler-combined"), "X axis"), choices = c("")),
                         selectInput(inputId = "ovparent", label = div(icon("sitemap"), "Parent"), choices = c("", "ungated")),
                         tags$style("#ovsamples {border-color:black}"),
                         actionButton(inputId = "ovsamples", label = "Select samples", icon("plus-square")),
                         br(), br(),
                         tags$style("#ovdisplayoptmain {border-color:black}"),
                         actionButton(inputId = "ovdisplayoptmain", label = div(icon("bars"), "Display options")),
                         width = 3),
                       mainPanel(
                         plotOutput("overlays"),
                         br(), br(), br(),
                         tags$style("#ovsampleorder {border-color:black}"),
                         actionButton(inputId = "ovsampleorder", label = "Edit order", icon("fas fa-pen")),
                         tags$style("#exportimageOverlay {border-color:black}"),
                         downloadButton(outputId = "exportimageOverlay", label = "Export image"),
                         width = 9),
              ),
              ##Proliferation----
              tabPanel(div(icon("chart-area"), "Proliferation"), value = "proliftab",
                       br(),
                       sidebarPanel(
                         tags$style("#step1 {border-color:black}"),
                         actionButton(inputId = "step1", label = div("Step 1: gate undivided", tags$i(class = "far fa-question-circle")) ),
                         br(), br(),
                         tags$style("#step2 {border-color:black}"),
                         actionButton(inputId = "step2", label = div("Step 2: apply model", tags$i(class = "far fa-question-circle")) ),
                         br(), br(),
                         tags$style("#step3 {border-color:black}"),
                         actionButton(inputId = "step3", label = div("Step 3: results", tags$i(class = "far fa-question-circle")) ),
                         br(), br(),
                         selectInput(inputId = "prolifsamplebutt", label = div(icon("vial"), "Sample"), choices = c("")),
                         checkboxInput(inputId = "proliflabel", label = "Numbers on top", value = TRUE),
                         checkboxInput(inputId = "prolifgrid", label = "Vertical lines", value = TRUE),
                         br(), br(),
                         tags$style("#applyprolif {border-color:black}"),
                         actionButton(inputId = "applyprolif", label = "Apply to Samples", icon("share-square")),
                         width = 3),
                       mainPanel(
                         htmlOutput("prolifstart"),
                         fluidRow(column(width = 3,
                                         div(style = "display:inline-block; position:fixed",
                                             actionButton(inputId = "prevprolifsample", width = "40px", icon("chevron-left"), style="color: white; background-color:black; border-color:black"),
                                             actionButton(inputId = "nextprolifsample", width = "40px", icon("chevron-right"), style="color: white; background-color:black; border-color:black") ) ),),
                         br(), br(),
                         fluidRow(column(width = 8,
                                         plotOutput("prolifplot", height = 450), ),
                                  column(width = 4,
                                         htmlOutput("proliftable"), ), ),
                         fluidRow(column(width = 8,
                                         tags$style("#exportimageProlif {border-color:black}"),
                                         downloadButton(outputId = "exportimageProlif", label = "Export image"), ),
                                  column(width = 4,
                                         tags$style("#exporttableProlif {border-color:black}"),
                                         downloadButton(outputId = "exporttableProlif", label = "Export complete table"), ), ),
                         width = 9)
              ),
              ##t-SNE----
              tabPanel(div(icon("chart-pie"), "t-SNE"), value = "tsnetab",
                       br(),
                       sidebarPanel(
                         selectInput(inputId = "tsnemode", label = div(icon("laptop-code"), "Mode"), choices = c("Heatmap", "Overlay Groups or Samples", "Overlay Populations")),
                         selectInput(inputId = "tsnehighlight", label = div(icon("highlighter"), "Highlight"), choices = c("")),
                         selectInput(inputId = "tsnegrouporsample", label = div(icon("vial"), "Events to show"), choices = c("All", "Group", "Sample")),
                         selectInput(inputId = "tsnegrouporsampleID", label = div(icon("vial"), "Group/Sample"), choices = c("")),
                         actionButton(inputId = "tsnegrouporsampleIDs", label = "Groups/Samples", icon("plus-square")),
                         actionButton(inputId = "tsnepopulations", label = "Populations", icon("plus-square")),
                         tags$style(shiny::HTML("[for=tsnegroups]+span>.irs>.irs-single, [for=tsnegroups]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
                         sliderInput("tsnegroups", label = "Number of groups", min = 1, max = 6, step = 1, value = 1, ticks = FALSE),
                         actionButton(inputId = "tsnesamples", label = "Select samples", icon("plus-square")),
                         br(), br(),
                         tags$style(shiny::HTML("[for=tsnedotsize]+span>.irs>.irs-single, [for=tsnedotsize]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
                         sliderInput("tsnedotsize", label = "Dot size", min = 4, max = 16, step = 1, value = 8, ticks = FALSE),
                         selectInput(inputId = "tsneparent", label = div(icon("sitemap"), "Parent"), choices = c("", "ungated")),
                         actionButton(inputId = "tsneparameters", label = "Select parameters", icon("plus-square")),
                         br(), br(),
                         selectInput(inputId = "tsneevents", label = div(icon("ellipsis-v"), "Events per sample"), choices = c(" " = 0, "1k" = 1000,  "5k" = 5000, "10k" = 10000, "25k" = 25000, "50k" = 50000)),
                         br(), br(),
                         actionButton(inputId = "tsnegenerate", label = "Generate t-SNE", icon("play-circle")),
                         width = 3),
                       mainPanel(
                         plotOutput("tsneplot"),
                         br(), br(), br(),
                         tags$style("#savetsneplot {border-color:black}"),
                         div(style="display:inline-block", downloadButton(outputId = "savetsneplot", label = "Export image")),
                         div(style="display:inline-block;vertical-align:top; width: 20px;", shiny::HTML("<br>")),
                         div(style="display:inline-block", textOutput("showingtsneevents")),
                         width = 9)
              ),
              ##Results----
              tabPanel(div(icon("table"), "Results"), value = "resulttab",
                       br(),
                       sidebarPanel(
                         selectInput(inputId = "rsstat", label = div(icon("percentage"), "Statistic"), choices = c("", "Freq. of parent", "Freq. of...", "Freq. of total", "Median", "Count")),
                         selectInput(inputId = "rsparent", label = div(icon("sitemap"), "Parent"), choices = c("")),
                         selectInput(inputId = "rspop", label = div(tags$i(class = "fas fa-clone"), "Population"), choices = c("")),
                         selectInput(inputId = "rsparameter", label = div(icon("ruler-combined"), "Parameter"), choices = c("")),
                         tags$style("#addresult {border-color:black}"),
                         actionButton(inputId = "addresult", label = "Add to table", icon("plus-square")),
                         tags$style("#editresult {border-color:black}"),
                         width = 3 ),
                       mainPanel(
                         tags$head(tags$style(shiny::HTML("#results td, #results th {border-color:black}")) ),
                         tableOutput("results"),
                         actionButton(inputId = "editresult", label = "Edit", icon("fas fa-pen")),
                         tags$style("#exporttable {border-color:black}"),
                         downloadButton(outputId = "exporttable", label = "Export table"),
                         width = 9 ) )
  ) )

#SERVER----
server <<- function(input, output, session) {
  shinyjs::runjs("$('#drawgatename').attr('maxlength',18)")
  shinyjs::runjs("$('#drawgatename').bind('keypress', function(event) {
    var regex = new RegExp('^[/]+$');
    if (regex.test(String.fromCharCode(!event.charCode ? event.which : event.charCode))) {
      event.preventDefault() } })")
  shinyjs::hide("savematrix")
  shinyjs::hide("cancelmatrix")
  shinyjs::hide("prolifsamplebutt")
  shinyjs::hide("proliflabel")
  shinyjs::hide("prolifgrid")
  shinyjs::hide("proliftable")
  shinyjs::hide("prevprolifsample")
  shinyjs::hide("nextprolifsample")
  shinyjs::hide("applyprolif")
  shinyjs::hide("exportimageProlif")
  shinyjs::hide("exporttableProlif")
  shinyjs::hide("tsnemode")
  shinyjs::hide("tsnedotsize")
  shinyjs::hide("tsnehighlight")
  shinyjs::hide("tsnegrouporsample")
  shinyjs::hide("tsnegrouporsampleID")
  shinyjs::hide("tsnegrouporsampleIDs")
  shinyjs::hide("tsnepopulations")
  shinyjs::disable("ovsamples")
  shinyjs::disable("step2")
  shinyjs::disable("step3")
  shinyjs::disable("exporttable")
  shinyjs::disable("editresult")
  shinyjs::disable("exportimageAncestry")
  shinyjs::disable("exportimageOverlay")
  shinyjs::disable("savetsneplot")
  shinyjs::disable("savefile")
  shinyjs::js$disableTab("plottab")
  shinyjs::js$disableTab("comptab")
  shinyjs::js$disableTab("ancestrytab")
  shinyjs::js$disableTab("overlaytab")
  shinyjs::js$disableTab("proliftab")
  shinyjs::js$disableTab("tsnetab")
  shinyjs::js$disableTab("resulttab")
  reactparent = shiny::reactiveValues(data = "root")
  reactnamechange = shiny::reactiveValues(data = NULL)
  reactrectangle = shiny::reactiveValues(data = FALSE); reactinterval = shiny::reactiveValues(data = FALSE); reactpolygon = shiny::reactiveValues(data = FALSE); reactquadrant = shiny::reactiveValues(data = FALSE)
  reacthover = shiny::reactiveValues(data = 0)
  hovercoords = shiny::reactiveValues(data = c(NULL))
  polygoncoords = shiny::reactiveValues(data = data.frame(NULL))
  plotactivator = shiny::reactiveValues(data = 0)
  hierarchyactivator = shiny::reactiveValues(data = 0)
  bgactivator = shiny::reactiveValues(data = 0)
  reactgatefont = shiny::reactiveValues(data = 5)
  reactshowgatename = shiny::reactiveValues(data = TRUE)
  reactaxisfont = shiny::reactiveValues(data = 15)
  reactshowaxis = shiny::reactiveValues(data = TRUE)
  ovactivator = shiny::reactiveValues(data = 0)
  showovplot = shiny::reactiveValues(data = FALSE)
  reactovsamples = shiny::reactiveValues(data = c())
  samplealertmessage = shiny::reactiveValues(data = "")
  reactmodaltitle = shiny::reactiveValues(data = "")
  reactsampleorder = shiny::reactiveValues(data = NULL)
  reactcomptable = shiny::reactiveValues(data = NULL)
  compplotactivator = shiny::reactiveValues(data = 0)
  currenttable <<- NULL
  reactreadnonly = shiny::reactiveValues(data = TRUE)
  appliedmatrix <<- "Cytometer-defined"
  reactaxiscustom = shiny::reactiveValues(data = NULL)
  prolifready = shiny::reactiveValues(data = FALSE)
  reacttsnepar = shiny::reactiveValues(data = NULL)
  reacttsnesamples = shiny::reactiveValues(data = NULL)
  reactgroupnames = shiny::reactiveValues(data = NULL)
  reacttsneparent = shiny::reactiveValues(data = NULL)
  tsneplotactivator = shiny::reactiveValues(data = 0)
  entiretSNE <<- NULL
  reacttsneIDs = shiny::reactiveValues(data = NULL)
  reacttsnepops = shiny::reactiveValues(data = NULL)
  reactovaxisfont = shiny::reactiveValues(data = 15)
  reactovshowaxis = shiny::reactiveValues(data = TRUE)
  ovboolean <<- FALSE
  tsneboolean <<- FALSE
  betweenpeaks <<- NULL
  refpeaks <<- NULL
  adjustedlabels <<- NULL
  concatsamples <<- NULL
  tsnelistofgroups <<- NULL
  tsnelistofsamples <<- NULL
  loadingfile <<- FALSE
  reactdir = shiny::reactiveValues(data = NULL)
  reactBCytofile = shiny::reactiveValues(data = NULL)
  reactheight = shiny::reactiveValues(data = 200)

  #Plot----
  ##left----
  observeEvent(input$customizeaxisY, {
    reactaxiscustom$data = input$channelY })
  observeEvent(input$customizeaxisX, {
    reactaxiscustom$data = input$channelX })
  observeEvent(reactaxiscustom$data, {
    setwd(filepath)
    poppaths = flowWorkspace::gs_get_pop_paths(gs)
    popparents = flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
    popgates = list()
    for(i in seq(poppaths)[-1]){
      popgates[[i]] = flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]] }
    popgates = popgates[-1]
    gs = flowWorkspace::gs_clone(uncompgs)
    flowCore::compensate(gs, compdfs[appliedmatrix][[1]])
    for(i in seq(popgates)){
      flowWorkspace::gs_pop_add(gs, popgates[[i]], parent = popparents[i]) }
    flowWorkspace::recompute(gs)
    predf = flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[match(input$sample, filelist)], reactparent$data)[[1]])@exprs[,reactaxiscustom$data]
    if(length(predf) > 0){
      dataf = data.frame(predf)
      if(nrow(dataf) > 5000){
        set.seed(4)
        tempdf <<- dataf[sample.int(nrow(dataf), 5000),]
      } else{
        tempdf <<- dataf }
      shiny::showModal(shiny::modalDialog(
        tags$style(shiny::HTML("[for=posslider]+span>.irs>.irs-single, [for=posslider]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
        tags$style(shiny::HTML("[for=negslider]+span>.irs>.irs-single, [for=negslider]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
        tags$style(shiny::HTML("[for=widthslider]+span>.irs>.irs-single, [for=widthslider]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
        fluidRow(column(width = 9, align = "right", plotOutput("axisplot", height = 440)),
                 column(width = 3, align = "left", style = "margin-top: 350px;",  actionButton(inputId = "resetbutton", label = "Reset", style="color: black; background-color:white; border-color:black"))),
        fluidRow(column(width = 12, align = "center", sliderInput("posslider", label = "Positive Decades", min = 2, max = 7, step = 0.02, value = currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][1], ticks = FALSE),
                        sliderInput("negslider", label = "Negative Decades", min = 0, max = 1, step = 0.1, value = currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][2], ticks = FALSE),
                        sliderInput("widthslider", label = "Width Basis", min = 0, max = 3, step = 0.1, value = log10(abs(currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][3])), ticks = FALSE))),
        footer = list(actionButton(inputId = "saveaxis", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
        easyClose = FALSE,
        size = "l"))
    } else{
      shinyjs::alert(paste0("The sample ", "'", filelist[match(input$sample, filelist)], "'", " has 0 events in the selected parent. Please choose another sample or change the parent."))
      reactaxiscustom$data = NULL } })

  observeEvent(input$resetbutton, {
    updateSliderInput(inputId = "posslider", value = 4.42)
    updateSliderInput(inputId = "negslider", value = 0)
    updateSliderInput(inputId = "widthslider", value = 2) })

  output$axisplot = renderPlot({
    par(mar = c(4,6,1,1) + 0.1, lwd = 2)
    widthBasis <<- -10^input$widthslider
    if(input$widthslider < 1){
      widthBasis <<- format(round(widthBasis, 2), nsmall = 2)
    } else if(input$widthslider < 2){
      widthBasis <<- format(round(widthBasis, 1), nsmall = 1)
    } else if(input$widthslider <= 3){
      widthBasis <<- round(widthBasis) }
    trans = flowWorkspace::flowjo_biexp(channelRange=4096, maxValue=262144, pos=input$posslider, neg=input$negslider, widthBasis=as.numeric(widthBasis))
    if(class(tempdf) == "data.frame"){
      tempdf <<- tempdf[[1]] }
    tempdf = trans(tempdf)
    Xlim = c(0, 4100)
    if(length(tempdf) > 1){
      referencehist = hist(tempdf, breaks = 20, plot = FALSE)
      refcounts = referencehist$counts
      refdensity = referencehist$density
      multiplier = refcounts/refdensity
      densline = density(tempdf)
      densline$y = densline$y*multiplier[1]
      maxdens = max(densline[2]$y)
      plot(densline, xaxt = "n", yaxt = "n", ann = FALSE, xlim = Xlim, ylim = c(0, maxdens*1.05), xaxs = "i", yaxs = "i")
      polygon(densline, col = "grey", border = "black", lwd = 1)
      axisticks(2, "histo")
    } else{
      plot(0, xaxt = "n", yaxt = "n", ann = FALSE, xaxs = "i", yaxs = "i", cex = 0) }
    title(ylab = "Count", cex.lab = 1+15/10, line = 3, font.lab = 2)
    title(xlab = names(ch)[which(colnames(gs) == reactaxiscustom$data)], cex.lab = 1+15/10, line = 3, font.lab = 2)
    axisticks(1, "log", c(input$posslider, input$negslider, as.numeric(widthBasis)))
  }, height = 440, width = 470)

  observeEvent(input$saveaxis, {
    if(!identical(currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]], c(input$posslider, input$negslider, as.numeric(widthBasis)))){
      inversetrans = flowWorkspace::flowjo_biexp(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][1], neg=currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][2], widthBasis=currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]][3], inverse = TRUE)
      trans = flowWorkspace::flowjo_biexp(channelRange=4096, maxValue=262144, pos=input$posslider, neg=input$negslider, widthBasis=as.numeric(widthBasis))
      currentcustomaxis[[which(fluochannels == reactaxiscustom$data)]] <<- c(input$posslider, input$negslider, as.numeric(widthBasis))
      poppaths = flowWorkspace::gs_get_pop_paths(gs)
      popparents <<- flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
      popgates = list()
      for(i in seq(poppaths)[-1]){
        popgates[[i]] = flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]] }
      popgates = popgates[-1]
      gs <<- flowWorkspace::gs_clone(uncompgs)
      flowCore::compensate(gs, compdfs[appliedmatrix][[1]])
      for(i in seq(popgates)){
        flowWorkspace::gs_pop_add(gs, popgates[[i]], parent = popparents[i]) }
      flowWorkspace::recompute(gs)
      names(currentcustomaxis) = fluochannels
      for(i in fluochannels[fluochannels != reactaxiscustom$data]){
        gs <<- transform(gs, transformerList(i, flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3]))) }
      gs <<- transform(gs, transformerList(reactaxiscustom$data, flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=input$posslider, neg=input$negslider, widthBasis=as.numeric(widthBasis))))
      flowWorkspace::recompute(gs)
      for(i in seq(poppaths)[-1]){
        gatechannels = names(flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]]@parameters)
        for(j in seq(gatechannels)){
          if(names(flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]]@parameters)[j] == reactaxiscustom$data){
            popgate = flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]]
            if(class(popgate)[1] == "polygonGate"){
              popgate@boundaries[,reactaxiscustom$data] = trans(inversetrans(popgate@boundaries[,reactaxiscustom$data]))
            } else{
              if(popgate@min[[j]] != "-Inf"){
                popgate@min[[j]] = trans(inversetrans(popgate@min[[j]]))
              } else{
                popgate@min[[j]] = -Inf }
              if(popgate@max[[j]] != "Inf"){
                popgate@max[[j]] = trans(inversetrans(popgate@max[[j]]))
              } else{
                popgate@max[[j]] = Inf } }
            popgatelist = list()
            for(k in seq(filelist)){
              popgatelist[[k]] = popgate }
            names(popgatelist) = filelist
            flowWorkspace::gs_pop_set_gate(gs, popgate@filterId, popgatelist)
            flowWorkspace::recompute(gs, popgate@filterId) } } }
      reactaxiscustom$data = NULL
      plotactivator$data = plotactivator$data + 1
    } else{
      reactaxiscustom$data = NULL }
    compplotactivator$data = compplotactivator$data + 1
    hierarchyactivator$data = hierarchyactivator$data + 1
    shiny::removeModal() })

  observeEvent(input$displayoptmain, {
    shiny::showModal(shiny::modalDialog(
      tags$style(shiny::HTML("[for=displayaxisfont]+span>.irs>.irs-single, [for=displayaxisfont]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
      checkboxInput("displayaxis", label = "Show axis titles", value = reactshowaxis$data),
      sliderInput("displayaxisfont", label = NULL, min = 10, max = 25, value = reactaxisfont$data, ticks = FALSE),
      br(),
      tags$style(shiny::HTML("[for=displaygatefont]+span>.irs>.irs-single, [for=displaygatefont]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
      checkboxInput("displaygatename", label = "Show gate name", value = reactshowgatename$data),
      sliderInput("displaygatefont", label = NULL, min = 1, max = 30, value = reactgatefont$data, ticks = FALSE),
      footer = list(actionButton(inputId = "reloadplotmodal", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color: white; border-color: black")),
      easyClose = TRUE,
      size = "s"))
    if(length(flowWorkspace::gs_get_pop_paths(gs)) > 1){
      shinyjs::show("displaygatename")
      shinyjs::show("displaygatefont")
    } else{
      shinyjs::hide("displaygatename")
      shinyjs::hide("displaygatefont") } })

  observeEvent(input$reloadplotmodal, {
    plotactivator$data = plotactivator$data + 1
    shiny::removeModal()
    reactshowaxis$data = input$displayaxis
    reactaxisfont$data = input$displayaxisfont
    reactshowgatename$data = input$displaygatename
    reactgatefont$data = input$displaygatefont })

  observeEvent(c(input$displayaxis, input$displayoptmain) , {
    if(length(input$displayaxis) != 0){
      if(input$displayaxis == TRUE){
        shinyjs::show("displayaxisfont")
      } else{
        shinyjs::hide("displayaxisfont") } } })

  observeEvent(input$help, {
    hierarchyactivator$data = hierarchyactivator$data + 1
    shiny::showModal(shiny::modalDialog(
      title = "BCyto (version 1.0)",
      "Helpful considerations will go here.",
      footer = NULL,
      size = "l",
      easyClose = TRUE)) })

  observeEvent(input$about, {
    shiny::showModal(shiny::modalDialog(
      title = "BCyto (version 1.0)",
      "Information about the software will go here.",
      footer = NULL,
      size = "l",
      easyClose = TRUE)) })

  ##main----
  observeEvent(input$type, {
    if(input$type != "Histogram"){
      if(input$channelY == ""){
        shiny::updateSelectInput(inputId = "channelY", selected = Y)
      } else{
        plotactivator$data = plotactivator$data + 1 }
      shinyjs::disable("interval"); shinyjs::enable("rectangle"); shinyjs::enable("polygon"); shinyjs::enable("quadrant")
    } else{
      plotactivator$data = plotactivator$data + 1
      shinyjs::enable("interval"); shinyjs::disable("rectangle"); shinyjs::disable("polygon"); shinyjs::disable("quadrant") } })

  observeEvent(input$nextsample, {
    if(match(input$sample, filelist) + 1 <= length(filelist)){
      shiny::updateSelectInput(inputId = "sample", selected = filelist[[match(input$sample, filelist) + 1]]) } })
  observeEvent(input$previoussample, {
    if(match(input$sample, filelist[onlysampleIDs]) > 1){
      shiny::updateSelectInput(inputId = "sample", selected = filelist[onlysampleIDs][[match(input$sample, filelist[onlysampleIDs]) - 1]]) } })

  observeEvent(c(input$sample, input$nextsample, input$previoussample), {
    if(input$sample != ""){
      if(match(input$sample, filelist[onlysampleIDs]) >= length(filelist[onlysampleIDs])){
        shinyjs::disable("nextsample")
      } else{
        shinyjs::enable("nextsample") }
      if(match(input$sample, filelist[onlysampleIDs]) <= 1){
        shinyjs::disable("previoussample")
      } else{
        shinyjs::enable("previoussample") } } })

  observeEvent(c(input$sample, input$channelY, input$channelX, reactparent$data, input$type, input$dotplotgateCANCEL), {
    reactrectangle$data = FALSE; reactinterval$data = FALSE; reactpolygon$data = FALSE; reactquadrant$data = FALSE })

  observeEvent(input$dotplotgateOK, {
    "%notin%" = Negate("%in%")
    if(input$drawgatename %notin% flowWorkspace::gs_get_pop_paths(gs, path = 1)){
      if(reactrectangle$data == TRUE){
        gate = flowCore::rectangleGate(.gate = matrix(c(input$mainplot_brush$xmin, input$mainplot_brush$xmax, input$mainplot_brush$ymin, input$mainplot_brush$ymax), ncol = 2, dimnames = list(c("min", "max"), c(input$channelX, input$channelY)))) }
      if(reactpolygon$data == TRUE){
        mat = matrix(ncol = 2, nrow = nrow(polygoncoords$data))
        for(i in seq(nrow(polygoncoords$data))){
          mat[i,] = c(polygoncoords$data$x[i], polygoncoords$data$y[i]) }
        colnames(mat) = c(input$channelX, input$channelY)
        gate = flowCore::polygonGate(mat) }
      if(reactinterval$data == TRUE){
        gate = flowCore::rectangleGate(.gate = matrix(c(input$mainplot_brush$xmin, input$mainplot_brush$xmax), ncol = 1, dimnames = list(c("min", "max"), c(input$channelX)))) }
      flowWorkspace::gs_pop_add(gs, gate, parent = isolate(reactparent$data), name = input$drawgatename)
      flowWorkspace::recompute(gs)
      detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, "OK", isolate(reactshowgatename$data), isolate(reactgatefont$data))
      updateTextInput(inputId = "drawgatename", value = "")
      reactrectangle$data = FALSE
      reactinterval$data = FALSE
      reactpolygon$data = FALSE
    } else{
      shinyjs::alert("Please choose a different gate name.") } })

  output$plotui = renderUI({
    if(input$type == "Histogram"){
      if(reactinterval$data == TRUE){
        brushopts = shiny::brushOpts("mainplot_brush", fill = "lightgrey", stroke = "black", opacity = 0.4, delay = 10, direction = "x")
      } else{
        brushopts = NULL }
      hoveropts = NULL
    } else{
      if(reactrectangle$data == TRUE){
        brushopts = shiny::brushOpts("mainplot_brush", fill = "lightgrey", stroke = "black", opacity = 0.4, delay = 10, direction = "xy")
        hoveropts = NULL
      } else{
        brushopts = NULL
        if(reactpolygon$data == TRUE || reactquadrant$data == TRUE){
          hoveropts = shiny::hoverOpts(id = "mainplot_hover", delay = 100, delayType = "debounce", nullOutside = TRUE)
        } else{
          hoveropts = NULL } } }
    plotOutput("mainplot", height = 440, click = "mainplot_click", dblclick = "mainplot_dblclick", brush = brushopts, hover = hoveropts) })

  output$mainplot = renderPlot({
    par(mar = c(4,6,1,1) + 0.1, lwd = 2)
    input$dotplotgateOK
    input$dotplotgateCANCEL
    plotactivator$data
    shinyjs::disable("editgate")
    reactnamechange$data = NULL
    session$resetBrush("mainplot_brush")
    if(isolate(reactrectangle$data) == FALSE && isolate(reactinterval$data) == FALSE && isolate(reactpolygon$data) == FALSE && isolate(reactquadrant$data) == FALSE){
      shinyjs::hide("dotplotgateOK", TRUE, "fade"); shinyjs::hide("dotplotgateCANCEL", TRUE, "fade"); shinyjs::hide("drawgatename", TRUE, "fade")
      shinyjs::delay(500, shinyjs::show("rectangle", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("polygon", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("quadrant", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("interval", TRUE, "fade"))
      updateTextInput(inputId = "drawgatename", value =  "", placeholder = "Please draw a gate")
      shinyjs::disable("drawgatename")
      polygoncoords$data = data.frame(NULL)
      if(isolate(input$type) == "Histogram"){
        shiny::updateSelectInput(inputId = "channelY", selected = "")
        shinyjs::disable("channelY")
      } else{
        shinyjs::enable("channelY")
        shiny::updateSelectInput(inputId = "channelY", selected = input$channelY) }
      newplot(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), isolate(reactshowaxis$data), isolate(reactaxisfont$data))
      detectgate(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data)) }
    if(isolate(reactrectangle$data) == TRUE){
      isolate({
        newplot(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), isolate(reactshowaxis$data), isolate(reactaxisfont$data))
        detectgate(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data)) }) }
    if(reactpolygon$data == TRUE){
      polygoncoords$data
      isolate({
        newplot(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), isolate(reactshowaxis$data), isolate(reactaxisfont$data))
        detectgate(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data)) })
      if(!is.null(isolate({input$mainplot_click}))){
        if(abs(isolate({input$mainplot_click$x})) > 0 && abs(isolate({input$mainplot_click$y})) > 0){
          for(i in seq(nrow(polygoncoords$data))){
            if(i > 1){
              segments(polygoncoords$data$x[i-1], polygoncoords$data$y[i-1], polygoncoords$data$x[i], polygoncoords$data$y[i], col = "red", lwd = 2) }
            points(polygoncoords$data$x[i], polygoncoords$data$y[i], pch = 19, cex = 2, col = "red") } } } }
    if(reactquadrant$data == TRUE){
      isolate({
        newplot(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), isolate(reactshowaxis$data), isolate(reactaxisfont$data))
        detectgate(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data)) })
      abline(v = hovercoords$data[1], h = hovercoords$data[2], lwd = 2, col = "red") }
    if(reactinterval$data == TRUE){
      isolate({
        newplot(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), isolate(reactshowaxis$data), isolate(reactaxisfont$data))
        detectgate(match(input$sample, filelist), input$channelX, input$channelY, isolate(reactparent$data), isolate(input$type), "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data)) }) }
    output$events = renderText({
      shownevents = flowWorkspace::gs_pop_get_stats(gs[input$sample])[,3][which(flowWorkspace::gs_get_pop_paths(gs) == isolate(reactparent$data))][[1]]
      totalevents = flowWorkspace::gs_pop_get_stats(gs[input$sample])[,3][1][[1]]
      paste0("Plotted events: ", format(shownevents, big.mark = ","), "/", format(totalevents, big.mark = ",")) })
    if(grepl("SC", input$channelX)){
      shinyjs::disable("customizeaxisX")
    } else{
      shinyjs::enable("customizeaxisX") }
    if(grepl("SC", input$channelY) || isolate(input$type) == "Histogram"){
      shinyjs::disable("customizeaxisY")
    } else{
      shinyjs::enable("customizeaxisY") }
  }, width = 470)

  observeEvent(c(input$rectangle, input$polygon, input$quadrant), {
    if(!is.null(gs)){
      detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, ">= 4", isolate(reactshowgatename$data), isolate(reactgatefont$data))
      if(length(detectedgate) >= 4){
        shiny::showModal(shiny::modalDialog(
          "Only 4 gates are allowed on a plot. Would you like to remove the existing gates from this plot?",
          footer = list(actionButton(inputId = "okmodal", label = "Delete existing gates", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
          easyClose = TRUE,
          size = "s"))
      } else{
        shinyjs::hide("rectangle", TRUE, "fade"); shinyjs::hide("polygon", TRUE, "fade"); shinyjs::hide("quadrant", TRUE, "fade"); shinyjs::hide("interval", TRUE, "fade") } } })

  observeEvent(input$rectangle, {
    detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, ">= 4", isolate(reactshowgatename$data), isolate(reactgatefont$data))
    if(length(detectedgate) < 4){
      reactrectangle$data = TRUE
      shinyjs::delay(500, shinyjs::show("dotplotgateOK", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("dotplotgateCANCEL", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("drawgatename", TRUE, "fade"))
      shinyjs::disable("dotplotgateOK") } })

  observeEvent(input$polygon, {
    detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, ">= 4", isolate(reactshowgatename$data), isolate(reactgatefont$data))
    if(length(detectedgate) < 4){
      reactpolygon$data = TRUE
      shinyjs::delay(500, shinyjs::show("dotplotgateOK", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("dotplotgateCANCEL", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("drawgatename", TRUE, "fade"))
      shinyjs::disable("dotplotgateOK") } })

  observeEvent(input$quadrant, {
    detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, ">= 4", isolate(reactshowgatename$data), isolate(reactgatefont$data))
    if(length(detectedgate) < 4){
      reactquadrant$data = TRUE
      shinyjs::delay(500, shinyjs::show("dotplotgateCANCEL", TRUE, "fade")) } })

  observeEvent(input$interval, {
    reactinterval$data = TRUE
    shinyjs::hide("rectangle", TRUE, "fade"); shinyjs::hide("polygon", TRUE, "fade"); shinyjs::hide("quadrant", TRUE, "fade"); shinyjs::hide("interval", TRUE, "fade")
    shinyjs::delay(500, shinyjs::show("dotplotgateOK", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("dotplotgateCANCEL", TRUE, "fade")); shinyjs::delay(500, shinyjs::show("drawgatename", TRUE, "fade"))
    shinyjs::disable("dotplotgateOK") })

  observeEvent(input$mainplot_hover, {
    if(!is.null(hovercoords$data[1])){
      if(max(input$mainplot_hover$x, hovercoords$data[1])/min(input$mainplot_hover$x, hovercoords$data[1]) > 1.1 || max(input$mainplot_hover$y, hovercoords$data[2])/min(input$mainplot_hover$y, hovercoords$data[2]) > 1.1){
        reacthover$data = reacthover$data + 1
        hovercoords$data[1] = input$mainplot_hover$x
        hovercoords$data[2] = input$mainplot_hover$y }
    } else{
      reacthover$data = reacthover$data + 1
      hovercoords$data[1] = input$mainplot_hover$x
      hovercoords$data[2] = input$mainplot_hover$y } })

  observeEvent(input$mainplot_click, {
    if(reactpolygon$data == TRUE){
      if(!is.null(polygoncoords$data)){
        if(nrow(polygoncoords$data) > 2){
          if(max(input$mainplot_click$x, polygoncoords$data$x[1])/min(input$mainplot_click$x, polygoncoords$data$x[1]) < 1.1 && max(input$mainplot_click$y, polygoncoords$data$y[1])/min(input$mainplot_click$y, polygoncoords$data$y[1]) < 1.1){
            polygoncoords$data = rbind(polygoncoords$data, c(polygoncoords$data$x[1], polygoncoords$data$y[1]))
            updateTextInput(inputId = "drawgatename", placeholder = "Type gate name")
            shinyjs::enable("drawgatename")
          } else{
            if(polygoncoords$data$x[nrow(polygoncoords$data)] != polygoncoords$data$x[1] && input$mainplot_click$y != polygoncoords$data$y[nrow(polygoncoords$data)]){
              polygoncoords$data = rbind(polygoncoords$data, c(input$mainplot_click$x, input$mainplot_click$y )) } }
        } else{
          polygoncoords$data = rbind(polygoncoords$data, c(input$mainplot_click$x, input$mainplot_click$y ))
          colnames(polygoncoords$data) = c("x", "y") } } }
    if(reactquadrant$data == TRUE){
      gate = flowCore::quadGate(.gate = matrix(c(input$mainplot_click$x, input$mainplot_click$y), ncol = 2, dimnames = list(c("value"), c(input$channelX, input$channelY))))
      flowWorkspace::gs_pop_add(gs, gate, parent = reactparent$data, names = c(paste0("Q1: ", input$channelX, "- ", input$channelY, "+"), paste0("Q2: ", input$channelX, "+ ", input$channelY, "+"), paste0("Q3: ", input$channelX, "+ ", input$channelY, "-"), paste0("Q4: ", input$channelX, "- ", input$channelY, "-")))
      flowWorkspace::recompute(gs)
      detectgate(1, input$channelX, input$channelY, reactparent$data, input$type, "OK", isolate(reactshowgatename$data), isolate(reactgatefont$data))
      updateTextInput(inputId = "drawgatename", value = "")
      hierarchyactivator$data = isolate(hierarchyactivator$data) + 1
      reactquadrant$data = FALSE } })

  observeEvent(c(input$mainplot_brush, input$mainplot_click, input$drawgatename), {
    if(reactrectangle$data == TRUE || reactinterval$data == TRUE){
      if(is.null(input$mainplot_brush)){
        shinyjs::disable("dotplotgateOK")
        updateTextInput(inputId = "drawgatename", value = "", placeholder = "Please draw a gate")
        shinyjs::disable("drawgatename")
      } else{
        if(input$drawgatename == ""){
          shinyjs::disable("dotplotgateOK")
          updateTextInput(inputId = "drawgatename", placeholder = "Type gate name")
          shinyjs::enable("drawgatename")
        } else{
          shinyjs::enable("dotplotgateOK") } } }
    if(reactpolygon$data == TRUE){
      if(input$drawgatename == ""){
        shinyjs::disable("dotplotgateOK")
      } else{
        shinyjs::enable("dotplotgateOK") } } })


  observeEvent(input$mainplot_dblclick, {
    allpops = flowWorkspace::gs_get_pop_paths(gs)[-1]
    detectedgateinfo = list()
    for(i in seq(allpops)){
      detectedgateinfo[[i]] = flowWorkspace::gs_pop_get_gate(gs, allpops[i])[1][[1]] }
    popparents = flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
    clickablegate = list()
    for(i in seq(detectedgateinfo)){
      channels = names(detectedgateinfo[[i]]@parameters)
      if(length(channels) == 1){
        channels[2] = Y }
      if(input$type != "Histogram"){
        if(channels[1] == input$channelX && channels[2] == input$channelY){
          clickablegate[[i]] = detectedgateinfo[[i]] }
      } else{
        if(channels[1] == input$channelX){
          clickablegate[[i]] = detectedgateinfo[[i]] } } }
    updateparent = FALSE
    for(i in clickablegate){
      if(class(i)[1] == "rectangleGate"){
        if(length(names(i@parameters)) > 1){
          if(input$mainplot_dblclick$x >= i@min[[1]] && input$mainplot_dblclick$x <= i@max[[1]] && input$mainplot_dblclick$y >= i@min[[2]] && input$mainplot_dblclick$y <= i@max[[2]]){
            updateparent = TRUE
            reactparent$data = poppaths[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == i@filterId)] }
        } else{
          if(input$mainplot_dblclick$x >= i@min[[1]] && input$mainplot_dblclick$x <= i@max[[1]]){
            updateparent = TRUE
            reactparent$data = poppaths[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == i@filterId)] } } }
      if(class(i)[1] == "polygonGate"){
        if(input$mainplot_dblclick$x >= min(i@boundaries[,1]) && input$mainplot_dblclick$x <= max(i@boundaries[,1]) && input$mainplot_dblclick$y >= min(i@boundaries[,2]) && input$mainplot_dblclick$y <= max(i@boundaries[,2])){
          updateparent = TRUE
          reactparent$data = poppaths[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == i@filterId)] } } }
    if(updateparent == TRUE){
      if(reactparent$data != "root"){
        channels = names(flowWorkspace::gs_pop_get_gate(gs, reactparent$data)[[1]]@parameters) }
      gateswiththisparent = detectedgateinfo[which(popparents == reactparent$data)]
      if(length(gateswiththisparent) > 0){
        channels = names(gateswiththisparent[[1]]@parameters) }
      shiny::updateSelectInput(inputId = "channelX", selected = channels[1])
      if(length(channels) == 2){
        if(input$type == "Histogram"){
          shiny::updateSelectInput(inputId = "type", selected = "Pseudocolor") }
        shiny::updateSelectInput(inputId = "channelY", selected = channels[2])
        if(input$channelX == channels[1] && input$channelY == channels[2]){
          plotactivator$data = plotactivator$data + 1 }
      } else{
        if(input$type != "Histogram"){
          shiny::updateSelectInput(inputId = "channelY", selected = "")
          shiny::updateSelectInput(inputId = "type", selected = "Histogram") }
        if(input$channelX == channels){
          plotactivator$data = plotactivator$data + 1 } } } })

  output$saveMainPlot = downloadHandler(
    filename = function(){
      paste0(substr(input$sample, 1, nchar(input$sample) - 4), ".png") },
    content = function(file){
      png(file, units = "in", height = 6, width = 6.44, res = 300)
      par(mar = c(4,6,1,1) + 0.1, lwd = 2)
      newplot(match(input$sample, filelist), input$channelX, input$channelY, reactparent$data, input$type, isolate(reactshowaxis$data), isolate(reactaxisfont$data))
      detectgate(match(input$sample, filelist), input$channelX, input$channelY, reactparent$data, input$type, "plot", isolate(reactshowgatename$data), isolate(reactgatefont$data))
      dev.off() } )

  ##right----
  output$hierarchy = renderPlot({
    input$dotplotgateOK
    hierarchyactivator$data
    input$okmodal
    if(length(flowWorkspace::gs_get_pop_paths(gs)) > 1){
      #BUG---- (cannot coerce type 'S4' to vector of type 'double')
      #other bug: the X is not set, and the Y is being FSC-A, and the axis is not being displayed on the X
      hplot = plot(gs)
      labels = flowWorkspace::gs_get_pop_paths(gs, path = 1)
      labels[1] = "ungated"
      for(i in seq(labels)){
        if(substr(labels[i], 1, 1) == "Q" && nchar(labels[i]) >= 13){
          labels[i] = substr(labels[i], 1, 2) } }
      names(labels) = nodes(hplot)
      nodeAttrs = list(label = labels)
      attrs = list(node = list(fillcolor = "white", shape = "box", width = 1, color = "gray90", style = "rounded"), graph = list(rankdir = "TB"))
      if(is.null(reactnamechange$data)){
        colour = "black"
        names(colour) = nodes(hplot)[which(flowWorkspace::gs_get_pop_paths(gs) == reactparent$data)]
      } else{
        if(reactparent$data == reactnamechange$data){
          colour = "red"
          names(colour) = nodes(hplot)[which(flowWorkspace::gs_get_pop_paths(gs) == reactparent$data)]
        } else{
          colour = c("black", "red")
          names(colour) = c(nodes(hplot)[which(flowWorkspace::gs_get_pop_paths(gs) == reactparent$data)], nodes(hplot)[which(flowWorkspace::gs_get_pop_paths(gs) == reactnamechange$data)]) } }
      nodeAttrs$color = colour
      hplot <<- plot(hplot, nodeAttrs = nodeAttrs, attrs = attrs); plot(hplot, nodeAttrs = nodeAttrs, attrs = attrs)
      shinyjs::js$enableTab("ancestrytab")
      shinyjs::js$enableTab("overlaytab")
      shinyjs::js$enableTab("proliftab")
      shinyjs::js$enableTab("tsnetab")
      shinyjs::js$enableTab("resulttab")
    } else{
      hplot <<- plot(agopen(new("graphNEL", nodes = c("ungated", "B"), edgemode = "directed", edgeL = list(ungated = "B", B = "ungated")), "", attrs = list(node = list(fillcolor = "white", shape = "rectangle", width = 1)), edgeAttrs = list(color = c("ungated~B" = "white")), nodeAttrs = list(color = c("B" = "white"), fontcolor = c("B" = "white"))))
      shinyjs::js$disableTab("ancestrytab")
      shinyjs::js$disableTab("overlaytab")
      shinyjs::js$disableTab("proliftab")
      shinyjs::js$disableTab("tsnetab")
      shinyjs::js$disableTab("resulttab") }
    pops = flowWorkspace::gs_get_pop_paths(gs, path = 1)[-1]
    if(loadingfile == FALSE){
      shiny::updateSelectInput(inputId = "bgpop", choices = c("", pops[-1]))
      shiny::updateSelectInput(inputId = "ovparent", choices = c("", pops)) }
    shiny::updateSelectInput(inputId = "tsneparent", choices = c("", pops))
    shiny::updateSelectInput(inputId = "rsparent", choices = c("", pops))
    loadingfile <<- FALSE })

  observeEvent(input$hierarchy_dblclick, {
    dataf = data.frame(x = c(1:length(hplot@AgNode)), y = c(1:length(hplot@AgNode)))
    for(i in seq(hplot@AgNode)){
      rownames(dataf)[i] = hplot@AgNode[[i]]@txtLabel@labelText
      dataf$x[i] = hplot@AgNode[[i]]@center@x
      dataf$y[i] = hplot@AgNode[[i]]@center@y }
    selected = shiny::nearPoints(dataf, input$hierarchy_dblclick, xvar = "x", yvar = "y", threshold = 25, maxpoints = 1, addDist = TRUE)
    if(length(rownames(selected)) > 0){
      allpops = flowWorkspace::gs_get_pop_paths(gs)[-1]
      detectedgateinfo = list()
      for(i in seq(allpops)){
        detectedgateinfo[[i]] = flowWorkspace::gs_pop_get_gate(gs, allpops[i])[1][[1]] }
      popparents = flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
      checkingpop = flowWorkspace::gs_get_pop_paths(gs)[which(rownames(dataf) == rownames(selected))]
      reactparent$data = checkingpop
      if(reactparent$data != "root"){
        channels = names(flowWorkspace::gs_pop_get_gate(gs, reactparent$data)[[1]]@parameters) }
      gateswiththisparent = detectedgateinfo[which(popparents == reactparent$data)]
      if(length(gateswiththisparent) > 0){
        channels = names(gateswiththisparent[[1]]@parameters) }
      shiny::updateSelectInput(inputId = "channelX", selected = channels[1])
      if(length(channels) == 2){
        if(input$type == "Histogram"){
          shiny::updateSelectInput(inputId = "type", selected = "Pseudocolor") }
        shiny::updateSelectInput(inputId = "channelY", selected = channels[2])
        if(input$channelX == channels[1] && input$channelY == channels[2]){
          plotactivator$data = plotactivator$data + 1 }
      } else{
        if(input$type != "Histogram"){
          shiny::updateSelectInput(inputId = "channelY", selected = "")
          shiny::updateSelectInput(inputId = "type", selected = "Histogram") }
        if(input$channelX == channels){
          plotactivator$data = plotactivator$data + 1 } } } })

  observeEvent(input$hierarchy_click, {
    dataf = data.frame(x = c(1:length(hplot@AgNode)), y = c(1:length(hplot@AgNode)))
    for(i in seq(hplot@AgNode)){
      rownames(dataf)[i] = hplot@AgNode[[i]]@txtLabel@labelText
      dataf$x[i] = hplot@AgNode[[i]]@center@x
      dataf$y[i] = hplot@AgNode[[i]]@center@y }
    selected = shiny::nearPoints(dataf, input$hierarchy_click, xvar = "x", yvar = "y", threshold = 25, maxpoints = 1, addDist = TRUE)
    if(length(rownames(selected)) > 0 && rownames(selected) != "ungated"){
      reactnamechange$data = flowWorkspace::gs_get_pop_paths(gs)[which(rownames(dataf) == rownames(selected))]
      shiny::updateSelectInput(inputId = "editgate", label = shiny::HTML("<span style='color: red'>Edit</span>"))
      shinyjs::enable("editgate")
    } else{
      shiny::updateSelectInput(inputId = "editgate", label = "Edit")
      shinyjs::disable("editgate")
      reactnamechange$data = NULL } })

  observeEvent(input$editgate, {
    shiny::showModal(shiny::modalDialog(
      paste0("Please type a new name for gate '", flowWorkspace::gs_get_pop_paths(gs, path = 1)[which(flowWorkspace::gs_get_pop_paths(gs) == reactnamechange$data)], "':"),
      br(), br(),
      textInput(inputId = "newnametextbox", label = NULL, width = "180px"),
      footer = list(actionButton(inputId = "changenamemodal", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black"), actionButton(inputId = "deletegate", label = "Delete gate", style="color: white; background-color:black; border-color:black")),
      easyClose = TRUE,
      size = "m")) })

  observeEvent(input$deletegate, {
    rootplotload = FALSE
    poppaths = flowWorkspace::gs_get_pop_paths(gs, path = 1)
    selectedgate = poppaths[which(flowWorkspace::gs_get_pop_paths(gs) == reactnamechange$data)]
    parent = poppaths[which(flowWorkspace::gs_get_pop_paths(gs) == reactparent$data)]
    if(selectedgate == parent){
      rootplotload = TRUE }
    for(i in strsplit(flowWorkspace::gs_get_pop_paths(gs), split = "/")[-1]){
      if(parent %in% i[-1]){
        if(selectedgate %in% i[-1]){
          if(which(i == selectedgate) < which(i == parent)){
            rootplotload = TRUE }
        } else{
          rootplotload = TRUE } } }
    if(rootplotload == TRUE){
      plotactivator$data = plotactivator$data + 1
      reactparent$data = "root"
    } else{
      plotactivator$data = plotactivator$data + 1 }
    flowWorkspace::gs_pop_remove(gs, selectedgate)
    hierarchyactivator$data = isolate(hierarchyactivator$data) + 1
    shiny::removeModal()
    shiny::updateSelectInput(inputId = "editgate", label = shiny::HTML("<span style='color: black'>Edit</span>")) })

  observeEvent(c(input$newnametextbox, input$editgate), {
    if(!is.null(input$newnametextbox)){
      if(input$newnametextbox == ""){
        shinyjs::disable("changenamemodal")
      } else{
        shinyjs::enable("changenamemodal") } } })

  observeEvent(input$changenamemodal, {
    if(input$newnametextbox %in% flowWorkspace::gs_get_pop_paths(gs, path = 1)){
      shinyjs::alert("Please choose a different gate name.")
    } else{
      flowWorkspace::gs_pop_set_name(gs, flowWorkspace::gs_get_pop_paths(gs, path = 1)[which(flowWorkspace::gs_get_pop_paths(gs) == reactnamechange$data)], input$newnametextbox)
      shiny::removeModal()
      plotactivator$data = plotactivator$data + 1 } })

  observeEvent(input$okmodal, {
    for(i in detectedgate){
      flowWorkspace::gs_pop_remove(gs, i) }
    shiny::removeModal()
    plotactivator$data = plotactivator$data + 1  })

  observeEvent(input$cancelmodal, {
    reactaxiscustom$data = NULL
    shiny::removeModal() })

  output$exportimageGates = downloadHandler(
    filename = "Gate hierarchy.png",
    content = function(file){
      png(file, units = "in", height = 6, width = 6.44, res = 300)
      par(mar = c(4,6,1,1) + 0.1, lwd = 2)
      plot(hplot)
      dev.off() } )

  #Compensation----
  ##left----
  observeEvent(input$creatematrix, {
    shiny::showModal(shiny::modalDialog(
      "Select the method.",
      br(), br(),
      if(is.null(compdfs["AutoSpill"][[1]])){
        actionButton(inputId = "createautospill", label = "AutoSpill*", style="color: white; background-color:black; border-color:black") },
      actionButton(inputId = "createmanually", label = "Manually", style="color: white; background-color:black; border-color:black"),
      br(), br(),
      if(is.null(compdfs["AutoSpill"][[1]])){
        "*Automated compensation calculation and matrix generation. This might take a few minutes" },
      footer = actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black"),
      easyClose = FALSE,
      size = "m")) })

  observeEvent(input$createautospill, {
    if(!is.null(compcontrolIDs) && length(compcontrolIDs) >= length(fluochannels) && all(sapply(filelist[compcontrolIDs], function(x) length(grep(x, list.files(filepath))) > 0))){
      shinyjs::hide("createmanually")
      shinyjs::hide("cancelmodal")
      shinyjs::disable("createautospill")
      shiny::withProgress(message = "Please wait...", detail = "", value = 0, max = 100, {
        asp = autospill::get.autospill.param()
        control.dir = getwd()
        control.def.file = paste0(getwd(), "/fcs_control.csv")
        onlycontrols = filelist[-sort(c(grep("Unstained", filelist), onlysampleIDs))]
        possiblechannelorder = list()
        cs = flowWorkspace::load_cytoset_from_fcs(filelist)
        for(i in seq(onlycontrols)){
          possibleparameter = flowWorkspace::cytoframe_to_flowFrame(cs[[onlycontrols[i]]])@description$`TUBE NAME`
          possibleparameter = strsplit(possibleparameter, " Stained Control")[[1]]
          possibleparameter = paste0(possibleparameter, "-A")
          possiblechannelorder[[i]] = strsplit(possibleparameter, " Stained Control")[[1]] }
        possiblechannelorder = unlist(possiblechannelorder)
        futurecsv = data.frame(possiblechannelorder, possiblechannelorder, possiblechannelorder, possiblechannelorder)
        colnames(futurecsv) = c("filename", "dye", "antigen", "wavelength")
        futurecsv[,1] = onlycontrols
        futurecsv[,3] = futurecsv[,4] = rep(NA, length(futurecsv[,3]))
        write.table(futurecsv, file = control.def.file, sep = ",", row.names = F)

        flow.control = suppressWarnings(autospill::read.flow.control(control.dir, control.def.file, asp))
        flow.gate = autospill::gate.flow.data(flow.control, asp)
        marker.spillover.unco.untr = autospill::get.marker.spillover(TRUE, flow.gate, flow.control, asp)
        refine.spillover.result = autospill::refine.spillover(marker.spillover.unco.untr, NULL, flow.gate, flow.control, asp) })
      rightorder = sapply(flow.control$marker.original, function(x) which(fluochannels == x))
      reorderedAutoSpilldf = cbind2(refine.spillover.result[[1]], rightorder)
      reorderedAutoSpilldf = reorderedAutoSpilldf[order(reorderedAutoSpilldf[,ncol(reorderedAutoSpilldf)]),]
      reorderedAutoSpilldf = reorderedAutoSpilldf[,-ncol(reorderedAutoSpilldf)]
      reorderedAutoSpilldf = rbind2(reorderedAutoSpilldf, rightorder)
      reorderedAutoSpilldf = reorderedAutoSpilldf[,order(reorderedAutoSpilldf[nrow(reorderedAutoSpilldf),])]
      reorderedAutoSpilldf = reorderedAutoSpilldf[-nrow(reorderedAutoSpilldf),]
      colnames(reorderedAutoSpilldf) = rownames(reorderedAutoSpilldf) = fluochannels
      compdfs["AutoSpill"][[1]] <<- reorderedAutoSpilldf
      shiny::removeModal()
      choices = names(compdfs)
      names(choices) = choices
      allnames = c()
      for(i in seq(choices)){
        if(choices[i] == appliedmatrix){
          allnames[i] = paste(appliedmatrix, "(applied)")
        } else{
          allnames[i] = choices[i] } }
      names(choices) = allnames
      shiny::updateSelectInput(inputId = "previewmatrix", choices = choices, selected = "AutoSpill")
    } else{
      shinyjs::alert("AutoSpill requires all compensation control files to be loaded.") } })

  observeEvent(input$createmanually, {
    shinyjs::disable("creatematrix")
    shinyjs::disable("applymatrix")
    shinyjs::show("savematrix")
    shinyjs::show("cancelmatrix")
    reactreadnonly$data = FALSE
    shiny::removeModal() })

  observeEvent(input$cancelmatrix, {
    shinyjs::hide("savematrix")
    shinyjs::hide("cancelmatrix")
    shinyjs::enable("creatematrix")
    reactreadnonly$data = TRUE
    compplotactivator$data = compplotactivator$data + 1 })

  observeEvent(input$savematrix, {
    shinyjs::hide("savematrix")
    shinyjs::hide("cancelmatrix")
    shinyjs::enable("creatematrix")
    reactreadnonly$data = TRUE
    compdfs[paste("Manual Matrix", length(grep("Manual Matrix", names(compdfs))) + 1)][[1]] <<- currenttable/100
    choices = names(compdfs)
    names(choices) = choices
    allnames = c()
    for(i in seq(choices)){
      if(choices[i] == appliedmatrix){
        allnames[i] = paste(appliedmatrix, "(applied)")
      } else{
        allnames[i] = choices[i] } }
    names(choices) = allnames
    shiny::updateSelectInput(inputId = "previewmatrix", choices = choices, selected = names(compdfs)[length(compdfs)]) })

  observeEvent(input$applymatrix,{
    shiny::showModal(shiny::modalDialog(
      paste0("The matrix '", input$previewmatrix, "' will be applied to all samples. Do you want to proceed?"),
      footer = list(actionButton(inputId = "OKapplymatrix", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "s")) })

  observeEvent(input$OKapplymatrix, {
    plotactivator$data = plotactivator$data + 1
    appliedmatrix <<- input$previewmatrix
    poppaths = flowWorkspace::gs_get_pop_paths(gs)
    popparents = flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
    popgates = list()
    for(i in seq(poppaths)[-1]){
      popgates[[i]] = flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]] }
    popgates = popgates[-1]
    gs <<- flowWorkspace::gs_clone(uncompgs)
    flowCore::compensate(gs, compdfs[appliedmatrix][[1]])
    for(i in seq(popgates)){
      flowWorkspace::gs_pop_add(gs, popgates[[i]], parent = popparents[i]) }
    flowWorkspace::recompute(gs)
    names(currentcustomaxis) = fluochannels
    for(i in seq(fluochannels)){
      gs = transform(gs, transformerList(fluochannels[i], flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3]))) }
    flowWorkspace::recompute(gs)
    choices = names(compdfs)
    names(choices) = choices
    allnames = c()
    for(i in seq(choices)){
      if(choices[i] == appliedmatrix){
        allnames[i] = paste(appliedmatrix, "(applied)")
      } else{
        allnames[i] = choices[i] } }
    names(choices) = allnames
    shiny::updateSelectInput(inputId = "previewmatrix", choices = choices, selected = appliedmatrix)
    shiny::removeModal()
    plotactivator$data })

  observeEvent(c(input$previewmatrix, input$cancelmatrix),{
    if(appliedmatrix == input$previewmatrix){
      shinyjs::disable("applymatrix")
    } else{
      shinyjs::enable("applymatrix") } })

  ##main----
  observeEvent(input$comp, {
    handsonobject = isolate(input$comp)
    invertedtable = unlist(handsonobject$data)
    table = vector(mode = "list", length = length(fluochannels))
    index = 0
    for(i in seq(fluochannels)){
      for(j in seq(fluochannels)){
        index = index + 1
        table[[j]][i] = invertedtable[[index]] } }
    table = as.data.frame(table, row.names = fluochannels)
    colnames(table) = rownames(table)
    if(!is.null(currenttable) && !identical(currenttable, table) && reactreadnonly$data == FALSE){
      compplotactivator$data = compplotactivator$data + 1 }
    currenttable <<- table })

  observe({
    output$compgraphs = renderPlot({
      input$previewmatrix
      compplotactivator$data
      if(isolate(reactreadnonly$data) == TRUE){
        compgraphs(match(input$previewsample, filelist), compdfs[input$previewmatrix][[1]], input$showuncomp)
      } else{
        compgraphs(match(input$previewsample, filelist), currenttable/100, input$showuncomp, compdfs[input$previewmatrix][[1]]) }
    }, height = reactheight$data, width = reactheight$data) })

  #Ancestry----
  ##left----
  observeEvent(input$bgtype, {
    if(isolate(input$bgtype) == "Backgating" && length(flowWorkspace::gs_get_pop_paths(gs)) > 2){
      shinyjs::enable("bgpop")
    } else{
      shiny::updateSelectInput(inputId = "bgpop", selected = "")
      shinyjs::disable("bgpop") } })

  observeEvent(input$bgpop, {
    if(input$bgtype == "Backgating"){
      bgactivator$data = isolate(bgactivator$data) + 1 } })

  ##main----
  output$ancestry = renderPlot({
    bgactivator$data
    ancestrygenerator <<- function(){
      par(mfrow = c(3, 5), mar = c(4,5,2,1) + 0.1, lwd = 2)
      if(input$bgtype != "" && isolate(input$bgsample) != ""){
        if(isolate(input$bgtype) == "Backgating" && isolate(input$bgpop) != ""){
          bgdf = as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[match(input$bgsample, filelist)], isolate(input$bgpop))[[1]])@exprs)
        } else{
          bgdf = NULL }
        gatechannels = list()
        for(i in seq(flowWorkspace::gs_get_pop_paths(gs))[-1]){
          gatechannels[[i-1]] = names(flowWorkspace::gs_pop_get_gate(gs[1], flowWorkspace::gs_get_pop_paths(gs)[i])[[1]]@parameters) }
        duplicateIDs = sort(unique(c(seq(gatechannels)[duplicated(gatechannels)], seq(gatechannels)[duplicated(gatechannels, fromLast = TRUE)])))
        duplicateparents = flowWorkspace::gs_pop_get_count_fast(gs[[1]], path = "auto")[,"Parent"][[1]][duplicateIDs]
        if(length(duplicateIDs) > 0){
          plotsequence = seq(gatechannels)[-duplicateIDs[duplicated(duplicateparents)]]
        } else{
          plotsequence = seq(gatechannels) }
        dat = data.frame(x = numeric(0), y = numeric(0))
        shiny::withProgress(message = "Generating", detail = "plot 0", value = 0, max = length(plotsequence)/10, {
          for(i in (plotsequence+1)){
            dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
            shiny::incProgress(0.1, detail = paste("plot", which((plotsequence+1) == i) ))
            Sys.sleep(0.1)
            gateinfo = flowWorkspace::gs_pop_get_gate(gs, flowWorkspace::gs_get_pop_paths(gs)[i])[[1]]
            Xch = names(gateinfo@parameters[1])
            if(length(gateinfo@parameters) > 1){
              Ych = names(gateinfo@parameters[2])
              type = input$bgtype
            } else{
              Ych = "Histogram"
              type = "Histogram" }
            gatefullname = flowWorkspace::gs_get_pop_paths(gs)[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == gateinfo@filterId)]
            if(i > 2){
              parent = flowWorkspace::gs_get_pop_paths(gs)[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == tail(strsplit(gatefullname, split = "/")[[1]], 2)[1])]
            } else{
              parent = "root" }
            newplot(match(input$bgsample, filelist), Xch, Ych, parent, type, TRUE, 2, bgdf, TRUE)
            detectgate(match(input$bgsample, filelist), Xch, Ych, parent, type, "plot", TRUE, 2, TRUE)
            if(parent == "root"){
              parent = "ungated"
            } else{
              parent = flowWorkspace::gs_get_pop_paths(gs, path = 1)[which(flowWorkspace::gs_get_pop_paths(gs) == parent)] }
            if(isolate(input$bgtype) == "Backgating" && isolate(input$bgpop) != ""){
              if(!is.null(bgdf)){
                if(parent == isolate(input$bgpop)){
                  title(main = list(parent, cex = 1.5, col = "red"))
                } else{
                  title(main = list(parent, cex = 1.5, col = "black")) } }
            } else{
              title(main = list(parent, cex = 1.5, col = "black")) } } })
        shinyjs::enable("exportimageAncestry")
      } else{
        shinyjs::disable("bgpop") } }
    ancestrygenerator()
  }, height = 480, width = 800)

  output$exportimageAncestry = downloadHandler(
    filename = "Ancestry plots.png",
    content = function(file){
      png(file, units = "in", height = 6.6, width = 11.05, res = 300)
      ancestrygenerator()
      dev.off() } )

  #Overlays----
  ##left----
  observeEvent(c(input$ovtype, input$ovtone, input$ovchannelY, input$ovchannelX, input$ovparent), {
    if(input$ovtype == "Overlaid histogram" || input$ovtype == "Offset histogram"){
      if(input$ovchannelY != ""){
        shiny::updateSelectInput(inputId = "ovchannelY", selected = "")
      } else{
        ovactivator$data = ovactivator$data + 1 }
      shinyjs::disable("ovchannelY")
      if(input$ovtype == "Overlaid histogram"){
        if(input$ovchannelX != "" && input$ovparent != ""){
          shinyjs::enable("ovsamples")
          if(length(isolate(reactovsamples$data)) > 2){
            reactovsamples$data = isolate(reactovsamples$data)[1:3] }
          if(length(isolate(reactovsamples$data)) < 3){
            showovplot$data = TRUE
          } else{
            showovplot$data = FALSE }
          showovplot$data = TRUE
        } else{
          shinyjs::disable("ovsamples")
          showovplot$data = FALSE }
        if(input$ovtone != ""){
          shiny::updateSelectInput(inputId = "ovtone", selected = "")
        } else{
          ovactivator$data = ovactivator$data + 1 }
        shinyjs::disable("ovtone")
      } else{
        if(input$ovtone != "" && input$ovchannelX != "" && input$ovparent != ""){
          shinyjs::enable("ovsamples")
          ovactivator$data = ovactivator$data + 1
          showovplot$data = TRUE
        } else{
          shinyjs::disable("ovsamples")
          showovplot$data = FALSE }
        shinyjs::enable("ovtone")
        ovactivator$data = ovactivator$data + 1 }
    } else{
      if(length(isolate(reactovsamples$data)) > 2){
        reactovsamples$data = isolate(reactovsamples$data)[1:2] }
      if(input$ovtype != "" && input$ovtone != "" && input$ovchannelY != "" && input$ovchannelX != "" && input$ovparent != ""){
        shinyjs::enable("ovsamples")
        ovactivator$data = ovactivator$data + 1
        showovplot$data = TRUE
      } else{
        shinyjs::disable("ovsamples")
        showovplot$data = FALSE }
      shinyjs::enable("ovtone")
      shinyjs::enable("ovchannelY")
      ovactivator$data = ovactivator$data + 1 } })

  observeEvent(input$ovsamples,{
    samplealertmessage$data = ""
    updateCheckboxGroupInput(inputId = "samplecheckbox", selected = isolate(reactovsamples$data))
    if(input$ovtype == "Dot plot"){
      modaltitle = "Select 2 samples for the Dot plot."
    } else if(input$ovtype == "Overlaid histogram"){
      modaltitle = "Select 2-3 samples for the Overlaid histogram."
    } else{
      modaltitle = "Select 2-10 samples for the Offset histogram." }
    shiny::showModal(shiny::modalDialog(
      tags$style("#selectedsamplealert{color: red}"),
      strong(modaltitle),
      checkboxGroupInput("samplecheckbox", label = "", choices = filelist[onlysampleIDs], width = "100%" ),
      footer = fluidRow(column(width = 8, align = "left", textOutput("selectedsamplealert")),
                        column(width = 4, align = "right", list(actionButton(inputId = "ovokbutton", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")), )),
      easyClose = TRUE,
      size = "m"))
    shinyjs::disable("ovokbutton") })

  observeEvent(input$samplecheckbox, {
    if(input$ovtype == "Dot plot"){
      if(length(input$samplecheckbox) == 2){
        samplealertmessage$data = ""
        shinyjs::enable("ovokbutton")
      } else{
        if(length(input$samplecheckbox) > 2){
          samplealertmessage$data = "Only 2 samples are allowed for Dot plots."
        } else{
          samplealertmessage$data = "" }
        shinyjs::disable("ovokbutton") }
    } else if(input$ovtype == "Overlaid histogram"){
      samplealertmessage$data = "Up to 3 samples are allowed for Overlaid histograms."
      if(length(input$samplecheckbox) == 2 || length(input$samplecheckbox) == 3){
        samplealertmessage$data = ""
        shinyjs::enable("ovokbutton")
      } else{
        if(length(input$samplecheckbox) > 3){
          samplealertmessage$data = "Up to 3 samples are allowed for Overlaid histograms."
        } else{
          samplealertmessage$data = "" }
        shinyjs::disable("ovokbutton") }
    } else{
      samplealertmessage$data = "Up to 10 samples are allowed for Overlaid histograms."
      if(length(input$samplecheckbox) > 1 && length(input$samplecheckbox) < 11){
        samplealertmessage$data = ""
        shinyjs::enable("ovokbutton")
      } else{
        if(length(input$samplecheckbox) > 10){
          samplealertmessage$data = "Up to 10 samples are allowed for Overlaid histograms."
        } else{
          samplealertmessage$data = "" }
        shinyjs::disable("ovokbutton") } } })

  output$selectedsamplealert = renderText({ samplealertmessage$data })

  observeEvent(input$ovokbutton, {
    ovactivator$data = ovactivator$data + 1
    reactovsamples$data = input$samplecheckbox
    shiny::removeModal() })

  observeEvent(input$ovdisplayoptmain, {
    shiny::showModal(shiny::modalDialog(
      tags$style(shiny::HTML("[for=ovdisplayaxisfont]+span>.irs>.irs-single, [for=ovdisplayaxisfont]+span>.irs-bar {background: black; border-top: black; border-bottom: black}")),
      checkboxInput("ovdisplayaxis", label = "Show axis titles", value = reactovshowaxis$data),
      sliderInput("ovdisplayaxisfont", label = NULL, min = 10, max = 25, value = reactovaxisfont$data, ticks = FALSE),
      footer = list(actionButton(inputId = "ovreloadplotmodal", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color: white; border-color: black")),
      easyClose = TRUE,
      size = "s")) })

  observeEvent(input$ovreloadplotmodal, {
    ovactivator$data = ovactivator$data + 1
    shiny::removeModal()
    reactovshowaxis$data = input$ovdisplayaxis
    reactovaxisfont$data = input$ovdisplayaxisfont })

  observeEvent(c(input$ovdisplayaxis, input$ovdisplayoptmain) , {
    if(length(input$ovdisplayaxis) != 0){
      if(input$ovdisplayaxis == TRUE){
        shinyjs::show("ovdisplayaxisfont")
      } else{
        shinyjs::hide("ovdisplayaxisfont") } } })

  ##main----
  output$overlays = renderPlot({
    ovactivator$data
    if(!is.null(reactovsamples$data)){
      if(isolate(showovplot$data) == TRUE){
        shinyjs::enable("ovdisplayoptmain")
        shinyjs::enable("ovsampleorder")
        shinyjs::enable("exportimageOverlay")
        if(isolate(input$ovparent) == "ungated"){
          currentparent = "root"
        } else{
          currentparent = isolate(input$ovparent) }
        ovboolean <<- TRUE
        overlay(match(isolate(reactovsamples$data), filelist), isolate(input$ovchannelX), isolate(input$ovchannelY), currentparent, isolate(input$ovtype), isolate(input$ovtone), isolate(reactovshowaxis$data), isolate(reactovaxisfont$data))
      } else{
        ovboolean <<- FALSE
        shinyjs::disable("ovdisplayoptmain")
        shinyjs::disable("ovsampleorder")
        shinyjs::disable("exportimageOverlay") }
    } else{
      ovboolean <<- FALSE
      shinyjs::disable("ovdisplayoptmain")
      shinyjs::disable("ovsampleorder")
      shinyjs::disable("exportimageOverlay") }
  }, height = 423, width = 800)

  observeEvent(input$ovsampleorder, {
    reactmodaltitle$data = 1
    reactsampleorder$data = NULL
    shiny::showModal(shiny::modalDialog(
      strong(textOutput("ordermodaltitle")),
      radioButtons("orderradiobutton", label = "", choices = rev(isolate(reactovsamples$data)), width = "100%", selected = "" ),
      footer = actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black"),
      easyClose = TRUE,
      size = "m")) })

  output$ordermodaltitle = renderText({ paste0("Select sample number ", reactmodaltitle$data, " (from top to bottom).") })

  observeEvent(input$orderradiobutton, {
    reactmodaltitle$data = reactmodaltitle$data + 1
    for(i in seq(input$orderradiobutton)){
      reactsampleorder$data = append(reactsampleorder$data, match(input$orderradiobutton, reactovsamples$data)) }
    if(length(reactsampleorder$data) != length(reactovsamples$data)){
      samplestoshow = reactovsamples$data
      for(i in rev(sort(reactsampleorder$data))){
        samplestoshow = samplestoshow[samplestoshow != samplestoshow[i]] }
      updateRadioButtons(inputId = "orderradiobutton", choices = rev(samplestoshow), selected = "")
    } else{
      reactovsamples$data = rev(reactovsamples$data[reactsampleorder$data])
      ovactivator$data = ovactivator$data + 1
      shiny::removeModal() } })

  output$exportimageOverlay = downloadHandler(
    filename = function(){
      paste0(input$ovtype, ".png") },
    content = function(file){
      png(file, units = "in", height = 5.8, width = 11.05, res = 300)
      if(isolate(input$ovparent) == "ungated"){
        currentparent = "root"
      } else{
        currentparent = isolate(input$ovparent) }
      overlay(match(isolate(reactovsamples$data), filelist), isolate(input$ovchannelX), isolate(input$ovchannelY), currentparent, isolate(input$ovtype), isolate(input$ovtone), isolate(reactovshowaxis$data), isolate(reactovaxisfont$data))
      dev.off() } )

  #Proliferation----
  ##left----
  observeEvent(input$applyprolif, {
    prolifready$data = TRUE })

  observeEvent(input$step1, {
    shiny::showModal(shiny::modalDialog(
      "Follow the instructions to start the proliferation tool.",
      footer = actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black"),
      easyClose = TRUE,
      size = "s")) })
  observeEvent(input$step2, {
    shiny::showModal(shiny::modalDialog(
      "Navigate through samples and choose the best fitting model to apply.",
      footer = actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black"),
      easyClose = TRUE,
      size = "s")) })
  observeEvent(input$step3, {
    shiny::showModal(shiny::modalDialog(
      "You can use the generated results or return to step 2 to apply a new model.",
      footer = list(actionButton(inputId = "gotostep2", label = "Apply a new model", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "s")) })
  observeEvent(input$gotostep2, {
    prolifready$data = FALSE
    shiny::removeModal() })

  ##main----
  output$prolifstart = renderUI({
    l0 = "This tool automates the detection and quantification of division peaks from cell proliferation assays."
    l1 = "To initialize the proliferation tool:"
    l2 = "-go to the Plot tab;"
    l3 = "-double click in the desired parent;"
    l4 = "-select the proliferation dye in the X axis (i.e. FITC-A :: CFSE);"
    l5 = "-select 'Histogram' in the Plot type;"
    l6 = "-make a narrow gate (that does not reach the X axis limit) on any sample for undivided cells (brightest peak);"
    l7 = "-certify that the population is named as 'undivided';"
    l8 = "-return to this tab."
    shiny::HTML("<br/>", paste(strong(l0)), "<br/>", "<br/>", paste(strong(l1),l2,l3,l4,l5,l6,l7,l8, sep = "<br/>")) })

  observeEvent(input$nextprolifsample, {
    if(match(input$prolifsamplebutt, filelist[onlysampleIDs]) + 1 <= length(filelist[onlysampleIDs])){
      shiny::updateSelectInput(inputId = "prolifsamplebutt", selected = filelist[onlysampleIDs][[match(input$prolifsamplebutt, filelist[onlysampleIDs]) + 1]]) } })
  observeEvent(input$prevprolifsample, {
    if(match(input$prolifsamplebutt, filelist[onlysampleIDs]) > 1){
      shiny::updateSelectInput(inputId = "prolifsamplebutt", selected = filelist[onlysampleIDs][[match(input$prolifsamplebutt, filelist[onlysampleIDs]) - 1]]) } })
  observeEvent(c(input$prolifsamplebutt, input$nextprolifsample, input$prevprolifsample), {
    if(input$prolifsamplebutt != ""){
      if(match(input$prolifsamplebutt, filelist[onlysampleIDs]) >= length(filelist[onlysampleIDs])){
        shinyjs::disable("nextprolifsample")
      } else{
        shinyjs::enable("nextprolifsample") }
      if(match(input$prolifsamplebutt, filelist[onlysampleIDs]) <= 1){
        shinyjs::disable("prevprolifsample")
      } else{
        shinyjs::enable("prevprolifsample") } } })

  output$prolifplot = renderPlot({
    input$tabs
    prolifgenerator <<- function(){
      undividedgate = flowWorkspace::gs_pop_get_gate(gs, "undivided")[[1]]
      prolifchannel <<- names(undividedgate@parameters)
      prolifparent <<- flowWorkspace::gs_pop_get_count_fast(gs[filelist[onlysampleIDs][1]], path = "auto")
      prolifparent <<- prolifparent[,3][which(prolifparent[,2] == "undivided")][[1]]
      par(mar = c(4,6,1,1) + 0.1, lwd = 2)
      newplot(match(input$prolifsamplebutt, filelist), prolifchannel, "Prolif", prolifparent, "Histogram", isolate(reactshowaxis$data), isolate(reactaxisfont$data))
      shinyjs::show("prolifsamplebutt")
      shinyjs::show("prevprolifsample")
      shinyjs::show("nextprolifsample")
      if(prolifready$data == FALSE){
        shinyjs::hide("proliflabel")
        shinyjs::hide("prolifgrid")
        shinyjs::hide("proliftable")
        shinyjs::hide("exportimageProlif")
        shinyjs::hide("exporttableProlif")
        shinyjs::disable("step1")
        shinyjs::disable("step3")
        shinyjs::show("prolifsamplebutt")
        shinyjs::show("prevprolifsample")
        shinyjs::show("nextprolifsample")
        shinyjs::show("applyprolif")
        shinyjs::enable("step2")
        undividedlimits = c(undividedgate@min[[1]], undividedgate@max[[1]])
        approxd0X = referencehist$breaks[referencehist$breaks > undividedlimits[1]]
        approxd0X = approxd0X[approxd0X < undividedlimits[2]]
        peakID = which(diff(sign(diff(densline$y))) == -2)
        Xpeaks = densline$x[peakID]
        Ypeaks = densline$y[peakID]
        betweenpeaks <<- densline$x[which(diff(sign(diff(densline$y))) == 2)]
        betweenpeaks <<- betweenpeaks[betweenpeaks < mean(approxd0X)*1.05]
        refpeaks <<- data.frame(x = NA, y = NA)
        for(i in seq(Xpeaks)){
          refpeaks[i,] <<- c(Xpeaks[[i]], Ypeaks[[i]]) }
        if(length(which(refpeaks[[1]] > mean(approxd0X)*1.05)) > 0){
          refpeaks <<- refpeaks[-which(refpeaks[[1]] > mean(approxd0X)*1.05),] }
        adjustedlabels <<- c()
        for(i in seq(nrow(refpeaks))){
          if(i == 1){
            adjustedlabels[i] <<- betweenpeaks[1] - mean(refpeaks[,1]-refpeaks[,1][1])/3
          } else if(i == nrow(refpeaks)){
            adjustedlabels[i] <<- refpeaks[,1][i]
          } else{
            adjustedlabels[i] <<- mean(c(betweenpeaks[i], betweenpeaks[i-1])) } }
        for(i in seq(refpeaks[,1])){
          text(rev(adjustedlabels)[i], maxdens*1.1, labels = i-1, cex = 1.5)
          abline(v = betweenpeaks[i], lty = 2) }
      } else{
        dfX = dataf[,prolifchannel]
        shinyjs::hide("applyprolif")
        shinyjs::disable("step1")
        shinyjs::disable("step2")
        shinyjs::show("exporttableProlif")
        shinyjs::show("exportimageProlif")
        shinyjs::show("proliflabel")
        shinyjs::show("prolifgrid")
        shinyjs::show("proliftable")
        shinyjs::enable("step3")
        reference = hist(dfX, breaks = 150, plot = FALSE)
        interval = mean(refpeaks[,1]-refpeaks[,1][1])/nrow(refpeaks)
        gatecoords <<- list()
        for(i in seq(refpeaks[,1])){
          if(input$proliflabel == FALSE){
            dX = c(refpeaks[,1][i]-interval, refpeaks[,1][i]+interval)
            dbreaks = reference$breaks[reference$breaks > dX[1]]
            dbreaks = dbreaks[dbreaks < dX[2]]
            dcounts = c()
            for(j in seq(dbreaks)){
              dcounts[j] = reference$counts[which(reference$breaks == dbreaks[j])] }
            y =  mean(dcounts)*13
            if(y > maxdens*1.1){
              y = maxdens*1.1 }
            text(adjustedlabels[i], y, labels = abs(i-nrow(refpeaks)), cex = 1.5)
          } else{
            text(rev(adjustedlabels)[i], maxdens*1.1, labels = i-1, cex = 1.5) }
          if(input$prolifgrid == TRUE){
            abline(v = betweenpeaks[i], lty = 2) }
          if(i == 1){
            gatecoords[[i]] <<- c(0, betweenpeaks[1])
          } else if(i == nrow(refpeaks)){
            gatecoords[[i]] <<- c(betweenpeaks[i-1], 4100)
          } else{
            gatecoords[[i]] <<- c(betweenpeaks[i-1], betweenpeaks[i]) } }
        divisionpercentages <<- c()
        divisioncounts <<- c()
        totalcells <<- length(dfX)
        for(i in seq(gatecoords)){
          filter = dfX[dfX > gatecoords[[i]][1]]
          filter = filter[filter < gatecoords[[i]][2]]
          divisioncounts[i] <<- length(filter)
          divisionpercentages[i] <<- length(filter)*100/totalcells }
        divisioncounts <<- rev(divisioncounts) #from d0 to d
        divisionpercentages <<- rev(as.numeric(format(round(divisionpercentages, 2), nsmall = 2))) #from d0 to d4
        totaldividedcount <<- totalcells-divisioncounts[1]
        totaldividedpercent <<- 100-divisionpercentages[1] } }
    poppaths = flowWorkspace::gs_get_pop_paths(gs, path = 1)
    if("undivided" %in% poppaths){
      prolifgenerator()
    } else{
      shinyjs::hide("prolifplot")
      shinyjs::hide("applyprolif")
      shinyjs::hide("exporttableProlif")
      shinyjs::hide("exportimageProlif")
      shinyjs::hide("proliflabel")
      shinyjs::hide("prolifgrid")
      shinyjs::hide("proliftable")
      shinyjs::disable("step2")
      shinyjs::disable("step3")
      shinyjs::hide("prolifsamplebutt")
      shinyjs::hide("prevprolifsample")
      shinyjs::hide("nextprolifsample")
      shinyjs::show("prolifstart")
      shinyjs::enable("step1") }
  }, height = 440, width = 470)

  output$exportimageProlif = downloadHandler(
    filename = function(){
      paste0(substr(input$prolifsamplebutt, 1, nchar(input$prolifsamplebutt) - 4), ".png") },
    content = function(file){
      png(file, units = "in", height = 6, width = 6.44, res = 300)
      poppaths = flowWorkspace::gs_get_pop_paths(gs, path = 1)
      if("undivided" %in% poppaths){
        prolifgenerator() }
      dev.off() } )

  ##right----
  output$proliftable = renderUI({
    input$tabs
    input$prolifsamplebutt
    proliftablegenerator = function(){
      l0 = "Proliferation assay data:"
      l1 = paste("Number of divisions:", (nrow(refpeaks) - 1))
      l2 = paste("Divided cells:", totaldividedcount)
      l3 = paste("Undivided cells:", divisioncounts[1])
      lines = c(l1, l2, l3)
      for(i in 2:nrow(refpeaks)){ # i = 2, 3, 4, 5
        lines[i+2] = paste0("Division ", (i - 1), ": ", divisioncounts[i]) }
      lines[length(lines) + 1] = paste("Percent divided:", as.numeric(format(round(totaldividedpercent, 2), nsmall = 2)) ) #this number needs to be converted into 2 decimals only
      lines[length(lines) + 1] = paste("Percent undivided:", divisionpercentages[1])
      constantlinesln = length(lines)
      for(i in 2:nrow(refpeaks)){ # i = 2, 3, 4, 5
        lines[i+(constantlinesln-1)] = paste0("Percent division ", (i - 1), ": ", divisionpercentages[i]) }
      shiny::HTML("<br/>", paste(strong(l0)), "<br/>", paste(lines, collapse = "<br/>", sep = "<br/>")) }
    if(prolifready$data == TRUE){
      proliftablegenerator()} })

  output$exporttableProlif = downloadHandler(
    filename = "Proliferation results.csv",
    content = function(file){
      columnnames = c("Number of divisions", "Divided cells", "Undivided cells")
      for(i in 2:nrow(refpeaks)){
        columnnames[i+2] = paste0("Division ", (i - 1)) }
      columnnames[length(columnnames) + 1] = "Percent divided"
      columnnames[length(columnnames) + 1] = "Percent undivided"
      constantcollenght = length(columnnames)
      for(i in 2:nrow(refpeaks)){
        columnnames[i+(constantcollenght-1)] = paste0("Percent division ", (i - 1)) }
      completetable = data.frame(matrix("", ncol = length(columnnames), nrow = length(filelist[onlysampleIDs])))
      for(j in seq(filelist[onlysampleIDs])){
        dfX = as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[filelist[onlysampleIDs][j]], prolifparent)[[1]])@exprs)[,prolifchannel]
        divisionpercentages = c()
        divisioncounts = c()
        totalcells = length(dfX)
        for(i in seq(gatecoords)){
          filter = dfX[dfX > gatecoords[[i]][1]]
          filter = filter[filter < gatecoords[[i]][2]]
          divisioncounts[i] = length(filter)
          divisionpercentages[i] = length(filter)*100/totalcells }
        divisioncounts = rev(divisioncounts)
        divisionpercentages = rev(as.numeric(format(round(divisionpercentages, 2), nsmall = 2)))
        totaldividedcount = totalcells-divisioncounts[1]
        totaldividedpercent = 100-divisionpercentages[1]
        columnvalues = c((nrow(refpeaks) - 1), totaldividedcount, divisioncounts[1])
        for(i in 2:nrow(refpeaks)){
          columnvalues[i+2] = divisioncounts[i] }
        columnvalues[length(columnvalues) + 1] = as.numeric(format(round(totaldividedpercent, 2), nsmall = 2))
        columnvalues[length(columnvalues) + 1] = divisionpercentages[1]
        constantcollenght = length(columnvalues)
        for(i in 2:nrow(refpeaks)){
          columnvalues[i+(constantcollenght-1)] = divisionpercentages[i] }
        completetable[j,] = columnvalues }
      colnames(completetable) = columnnames
      rownames(completetable) = filelist[onlysampleIDs]
      write.csv(completetable, file) } )

  #t-SNE----
  ##left----
  observeEvent(input$tsnegroups, {
    if(length(reacttsnesamples$data) != input$tsnegroups){
      reacttsnesamples$data = NULL } })

  observeEvent(input$tsnesamples,{
    samplecolum = list()
    for(i in seq(filelist[onlysampleIDs])){
      samplecolum[[i]] = list(br(), br(), filelist[onlysampleIDs][i]) }
    groupcolums = list()
    for(i in seq(input$tsnegroups)){
      groupcolums[[i]] = list() }
    if(!is.null(reacttsnesamples$data)){
      for(j in seq(input$tsnegroups)){
        for(i in seq(filelist[onlysampleIDs])){
          groupcolums[[j]][[i]] = checkboxInput(inputId = paste0("group", j, "sample", i), label = "", value = reacttsnesamples$data[[j]][[i]]) } }
    } else{
      for(j in seq(input$tsnegroups)){
        for(i in seq(filelist[onlysampleIDs])){
          groupcolums[[j]][[i]] = checkboxInput(inputId = paste0("group", j, "sample", i), label = "") } } }
    checkboxcolumns <<- list()
    for(i in seq(input$tsnegroups)){
      if(!is.null(reactgroupnames$data) && length(reactgroupnames$data) == input$tsnegroups){
        textgroup = textInput(paste0("textgroup", i), label = "", value = reactgroupnames$data[[i]])
      } else{
        textgroup = textInput(paste0("textgroup", i), label = "") }
      checkboxcolumns[[i]] <<- column(width = 1, align = "center", style = "margin-top: 10px;", list(strong(paste("Group", i)), textgroup, groupcolums[[i]])) }
    shiny::showModal(shiny::modalDialog(
      strong("Select samples by group for concatenation."),
      br(),
      fluidRow(column(width = 5, align = "right", list(br(), br(), br(), br(), samplecolum)),
               checkboxcolumns ),
      footer = list(actionButton(inputId = "tsneoksampbutton", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "l")) })

  observeEvent(input$tsneoksampbutton,{
    reacttsnesamples$data = NULL
    reactgroupnames$data = NULL
    for(i in seq(input$tsnegroups)){
      reactgroupnames$data[[i]] = input[[paste0("textgroup", i)]]
      reacttsnesamples$data[[i]] = list() }
    for(j in seq(input$tsnegroups)){
      for(i in seq(filelist[onlysampleIDs])){
        reacttsnesamples$data[[j]][[i]] = input[[paste0("group", j, "sample", i)]] } }
    concatsamples <<- list()
    for(j in seq(input$tsnegroups)){
      concatsamples[[j]] <<- list() }
    for(j in seq(input$tsnegroups)){
      concatsamples[[j]] <<- which(reacttsnesamples$data[[j]] == TRUE) }
    if(length(concatsamples[[1]]) == 0){
      reacttsnesamples$data = NULL
      reactgroupnames$data = NULL
    } else{
      for(i in seq(concatsamples)){
        concatsamples[[i]] <<- sapply(concatsamples[[i]], function(x) which(filelist == filelist[onlysampleIDs][x])) } }
    if(length(grep(TRUE, unlist(lapply(concatsamples, function(x) length(x) == 0)))) > 0){
      concatsamples <<- list()
      shinyjs::alert("Select at least one sample per group.")
    } else{
      if(is.null(reactgroupnames$data)){
        concatsamples <<- list()
        shinyjs::alert("Please type the name of each group in the corresponding text boxes.")
      } else{
        if(length(grep(TRUE, sapply(reactgroupnames$data, function(x) x != ""))) < length(concatsamples)){
          concatsamples <<- list()
          shinyjs::alert("Please type the name of each group in the corresponding text boxes.")
        } else if(length(reactgroupnames$data) > 1 && length(unique(reactgroupnames$data)) < length(reactgroupnames$data)){
          shinyjs::alert("Group names cannot be identical.")
        } else{
          shiny::removeModal() } } } })

  observeEvent(input$tsneparameters,{
    updateCheckboxGroupInput(inputId = "tsneparcheckbox", selected = isolate(reacttsnepar$data))
    shiny::showModal(shiny::modalDialog(
      strong("Select at least 2 parameters to include in the analysis"),
      checkboxGroupInput("tsneparcheckbox", label = "", choices = unlist(ch)[unlist(ch) != "Time"], width = "100%" ),
      footer = list(actionButton(inputId = "tsneokparbutton", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "m")) })

  observeEvent(input$tsneokparbutton,{
    if(length(input$tsneparcheckbox) < 2){
      shinyjs::alert("Please select at least 2 parameters.")
    } else{
      reacttsnepar$data = input$tsneparcheckbox
      shiny::removeModal() } })

  observeEvent(c(reacttsnesamples$data, input$tsneparent, reacttsnepar$data, input$tsneevents), {
    if(!is.null(reacttsnesamples$data) && input$tsneparent != "" && !is.null(reacttsnepar$data) && input$tsneevents > 0){
      shinyjs::enable("tsnegenerate")
    } else{
      shinyjs::disable("tsnegenerate") } })

  observeEvent(input$tsnegenerate, {
    if(is.null(entiretSNE)){
      if(all(sapply(unlist(concatsamples), function(x) nrow(flowWorkspace::gs_pop_get_data(gs[x], input$tsneparent))[[1]] >= as.numeric(input$tsneevents)))){
        minutes = as.numeric(format(round(length(unlist(concatsamples))*as.numeric(input$tsneevents)/1000*2.5/60, 1), nsmall = 1))
        shiny::showModal(shiny::modalDialog(
          "The t-SNE tool might require longer processing time than other tools in this software.",
          br(), br(),
          "Before proceeding, make sure that your CPU is not being used by heavy processes.",
          br(), br(),
          "Concatenation of the selected samples/number of events and t-SNE generation are estimated to take ", strong(minutes, "-", minutes*2, " minutes."),
          br(), br(),
          textOutput("question"),
          footer = list(actionButton(inputId = "tsnegenerateok", label = "Generate t-SNE", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
          easyClose = FALSE,
          size = "m"))
      } else{
        lowernumber = min(sapply(unlist(concatsamples), function(x) nrow(flowWorkspace::gs_pop_get_data(gs[x], input$tsneparent))[[1]]))
        if(lowernumber < 1000){
          shinyjs::alert("There are not enough events (<1K) in at least one of the selected samples/parent. Try to change the samples or the Parent.")
        } else{
          if(lowernumber < 5000){
            maxevents = "1K"
          } else if(lowernumber < 10000){
            maxevents = "5K"
          } else if(lowernumber < 25000){
            maxevents = "10K"
          } else if(lowernumber < 50000){
            maxevents = "25K" }
          shinyjs::alert(paste("At least one selected samples/parent don't have enough events. For the selected samples/parent, only", maxevents, "events are possible.")) } }
    } else{
      shiny::showModal(shiny::modalDialog(
        "Do you want to make a new t-SNE?",
        footer = list(actionButton(inputId = "tsnemakenew", label = "Make new t-SNE", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
        easyClose = TRUE,
        size = "m")) } })

  output$question = renderText("Do you want to proceed?")

  observeEvent(input$tsnemakenew,{
    entiretSNE <<- NULL
    shinyjs::hide("tsnemode")
    shinyjs::hide("tsnedotsize")
    shinyjs::hide("tsnehighlight")
    shinyjs::hide("tsnegrouporsample")
    shinyjs::hide("tsnegrouporsampleID")
    shinyjs::hide("savetsneplot")
    shinyjs::hide("tsnegrouporsampleIDs")
    shinyjs::hide("tsnepopulations")
    shinyjs::hide("showingtsneevents")
    shinyjs::disable("tsnegenerate")
    shinyjs::show("tsnegroups")
    shinyjs::show("tsnesamples")
    shinyjs::show("tsneparent")
    shinyjs::show("tsneparameters")
    shinyjs::show("tsneevents")
    tsneplotactivator$data = tsneplotactivator$data + 1
    shiny::removeModal() })

  ##main----
  observeEvent(input$tsnegenerateok, {
    setwd(filepath)
    reacttsneparent$data = isolate(input$tsneparent)
    populations = gh_pop_get_descendants(gs[[1]], isolate(input$tsneparent), path = 1)
    for(i in seq(length(unlist(concatsamples)))){
      if(i == 1){
        preconcat = as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[unlist(concatsamples)[i]], "root")[[1]])@exprs)
        for(k in c(1, seq(populations) + 1)){
          if(k < length(c(1, seq(populations) + 1))){
            preconcat = cbind2(preconcat, gh_pop_get_indices(gs[[unlist(concatsamples)[i]]], populations[k]))
          } else{
            preconcat = cbind2(preconcat, gh_pop_get_indices(gs[[unlist(concatsamples)[i]]], isolate(input$tsneparent))) } }
        preconcat = preconcat[which(preconcat[ncol(preconcat)][[1]] == TRUE),]
        preconcat = preconcat[,-length(preconcat)]
        preconcat = preconcat[sample.int(nrow(preconcat), as.numeric(isolate(input$tsneevents))),]
      } else{
        temppreconcat = as.data.frame.array(flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(gs[unlist(concatsamples)[i]], "root")[[1]])@exprs)
        for(k in c(1, seq(populations) + 1)){
          if(k < length(c(1, seq(populations) + 1))){
            temppreconcat = cbind2(temppreconcat, gh_pop_get_indices(gs[[unlist(concatsamples)[i]]], populations[k]))
          } else{
            temppreconcat = cbind2(temppreconcat, gh_pop_get_indices(gs[[unlist(concatsamples)[i]]], isolate(input$tsneparent))) } }
        temppreconcat = temppreconcat[which(temppreconcat[ncol(temppreconcat)][[1]] == TRUE),]
        temppreconcat = temppreconcat[,-length(temppreconcat)]
        temppreconcat = temppreconcat[sample.int(nrow(temppreconcat), as.numeric(isolate(input$tsneevents))),]
        preconcat = rbind2(preconcat, temppreconcat) } }
    for(i in seq(fluochannels)){
      trans = flowWorkspace::flowjo_biexp(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3], inverse = TRUE)
      preconcat[,fluochannels[i]] = trans(preconcat[,fluochannels[i]]) }
    concat <<- NULL
    for(i in seq(isolate(input$tsneparcheckbox))){
      if(i == 1){
        concat <<- preconcat[,isolate(input$tsneparcheckbox)[1]]
      } else{
        concat <<- cbind2(concat, preconcat[,isolate(input$tsneparcheckbox)[i]]) } }
    for(i in rev(ncol(preconcat) - (seq(populations) - 1))){
      concat <<- cbind2(concat, preconcat[,i]) }
    concat <<- as.data.frame(concat)
    colnames(concat) <<- c(isolate(input$tsneparcheckbox), populations)
    shinyjs::hide("tsnegenerateok")
    shinyjs::hide("cancelmodal")
    shinyjs::hide("question")
    shiny::withProgress(message = "Generating t-SNE...", detail = "", value = 0, max = 100, {
      entiretSNE <<- Rtsne::Rtsne(as.data.frame(concat[,1:length(isolate(input$tsneparcheckbox))]), check_duplicates = FALSE, verbose = TRUE)
      entiretSNE <<- entiretSNE$Y })
    shinyjs::hide("tsnegroups")
    shinyjs::hide("tsnesamples")
    shinyjs::hide("tsneparent")
    shinyjs::hide("tsneparameters")
    shinyjs::hide("tsneevents")
    shinyjs::show("tsnemode")
    shinyjs::show("tsnedotsize")
    shinyjs::show("tsnehighlight")
    shinyjs::show("tsnegrouporsample")
    shinyjs::show("tsnegrouporsampleID")
    shinyjs::show("showingtsneevents")
    shinyjs::show("savetsneplot")
    shiny::removeModal()
    tsnelistofgroups <<- list()
    for(i in seq(concatsamples)){
      tsnelistofgroups[[i]] <<- i
      names(tsnelistofgroups)[[i]] <<- reactgroupnames$data[[i]] }
    tsnelistofsamples <<- list()
    for(i in seq(length(unlist(concatsamples)))){
      tsnelistofsamples[[i]] <<- i
      names(tsnelistofsamples)[[i]] <<- filelist[unlist(concatsamples)[i]] }
    if(length(tsnelistofgroups) < 2 && length(tsnelistofsamples) < 2){
      shiny::updateSelectInput(inputId = "tsnemode", choices = c("Heatmap", "Overlay Populations"))
    } else{
      shiny::updateSelectInput(inputId = "tsnemode", choices = c("Heatmap", "Overlay Groups or Samples", "Overlay Populations")) }
    tempchannels = unlist(ch)[unlist(ch) != "Time"]
    availableparameters <<- list()
    for(i in seq(input$tsneparcheckbox)){
      availableparameters[[i]] <<- which(tempchannels == input$tsneparcheckbox[i]) }
    availableparameters <<- tempchannels[unlist(availableparameters)]
    shiny::updateSelectInput(inputId = "tsnehighlight", choices = c(availableparameters))
    shiny::updateSelectInput(inputId = "tsnemode", selected = "Heatmap")
    shiny::updateSelectInput(inputId = "tsnegrouporsample", selected = "All")
    shiny::updateSelectInput(inputId = "tsnegrouporsampleID", selected = "") })

  observeEvent(input$tsnegrouporsample, {
    if(exists("concatsamples")){
      if(input$tsnemode == "Heatmap" || input$tsnemode == "Overlay Populations"){
        if(input$tsnegrouporsample == "All"){
          shinyjs::disable("tsnegrouporsampleID")
          shiny::updateSelectInput(inputId = "tsnegrouporsampleID", choices = c(""), selected = "")
        } else{
          shinyjs::enable("tsnegrouporsampleID")
          if(input$tsnegrouporsample == "Group"){
            shiny::updateSelectInput(inputId = "tsnegrouporsampleID", choices = tsnelistofgroups, selected = 1)
            if(isolate(input$tsnegrouporsampleID) == 1){
              tsneplotactivator$data = tsneplotactivator$data + 1 }
          } else if(input$tsnegrouporsample == "Sample"){
            shiny::updateSelectInput(inputId = "tsnegrouporsampleID", choices = tsnelistofsamples, selected = 1)
            if(isolate(input$tsnegrouporsampleID) == 1){
              tsneplotactivator$data = tsneplotactivator$data + 1 } } }
      } else if(input$tsnemode == "Overlay Groups or Samples"){
        if(input$tsnegrouporsample != "All"){
          shinyjs::enable("tsnegrouporsampleIDs")
        } else{
          shinyjs::disable("tsnegrouporsampleIDs") }
        reacttsneIDs$data = NULL
        tsneplotactivator$data = tsneplotactivator$data + 1 } } })

  observeEvent(input$tsnemode,{
    if(!is.null(entiretSNE)){
      if(input$tsnemode == "Heatmap" || input$tsnemode == "Overlay Populations"){
        if(input$tsnegrouporsample != "All"){
          if(input$tsnegrouporsample == "Group"){
            choices = tsnelistofgroups
          } else if(input$tsnegrouporsample == "Sample"){
            choices = tsnelistofsamples }
          if(input$tsnegrouporsampleID != 1){
            shiny::updateSelectInput(inputId = "tsnegrouporsampleID", choices = choices, selected = 1)
          } else{
            tsneplotactivator$data = tsneplotactivator$data + 1 }
          shinyjs::enable("tsnegrouporsampleID")
        } else{
          if(input$tsnegrouporsampleID != ""){
            shiny::updateSelectInput(inputId = "tsnegrouporsampleID", choices = c(""), selected = "")
          } else{
            tsneplotactivator$data = tsneplotactivator$data + 1 }
          shinyjs::disable("tsnegrouporsampleID") }
      } else if(input$tsnemode == "Overlay Groups or Samples"){
        if(input$tsnegrouporsample != "All"){
          shinyjs::enable("tsnegrouporsampleIDs")
        } else{
          shinyjs::disable("tsnegrouporsampleIDs") } }
      if(input$tsnemode == "Heatmap"){
        shinyjs::hide("tsnepopulations")
        shinyjs::hide("tsnegrouporsampleIDs")
        shinyjs::show("tsnegrouporsampleID")
        shinyjs::show("tsnegrouporsample")
        shinyjs::show("tsnehighlight")
      } else if(input$tsnemode == "Overlay Groups or Samples"){
        shinyjs::hide("tsnehighlight")
        shinyjs::hide("tsnegrouporsampleID")
        shinyjs::hide("tsnepopulations")
        shinyjs::show("tsnegrouporsampleIDs")
        reacttsneIDs$data = NULL
        tsneplotactivator$data = tsneplotactivator$data + 1
      } else if(input$tsnemode == "Overlay Populations"){
        shinyjs::hide("tsnehighlight")
        shinyjs::hide("tsnegrouporsampleIDs")
        shinyjs::show("tsnegrouporsampleID")
        shinyjs::show("tsnepopulations") } } })

  observeEvent(input$tsnegrouporsampleIDs,{
    if(is.null(reacttsneIDs$data)){
      if(input$tsnegrouporsample == "Group"){
        listofgrouporsamples = tsnelistofgroups
      } else if(input$tsnegrouporsample == "Sample"){
        listofgrouporsamples = tsnelistofsamples }
    } else{
      if(input$tsnegrouporsample == "Group"){
        listofgrouporsamples = tsnelistofgroups
      } else if(input$tsnegrouporsample == "Sample"){
        listofgrouporsamples = tsnelistofsamples } }
    shiny::showModal(shiny::modalDialog(
      "Select 2 groups/samples",
      checkboxGroupInput("tsneIDs", label = "", choices = listofgrouporsamples, selected = reacttsneIDs$data),
      footer = list(actionButton(inputId = "tsnesettingIDs", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "m"))
    if(length(input$tsneIDs) == 0){
      shinyjs::disable("tsnesettingIDs")
    } else{
      if(length(input$tsneIDs) != 2){
        shinyjs::disable("tsnesettingIDs")
      } else{
        shinyjs::enable("tsnesettingIDs") } } })

  observeEvent(input$tsneIDs,{
    if(input$tsneIDs[1] == ""){
      shinyjs::disable("tsnesettingIDs")
    } else{
      if(length(input$tsneIDs) != 2){
        shinyjs::disable("tsnesettingIDs")
      } else{
        shinyjs::enable("tsnesettingIDs") } } })

  observeEvent(input$tsnesettingIDs,{
    reacttsneIDs$data = input$tsneIDs
    shiny::removeModal() })

  observeEvent(input$tsnepopulations,{
    popoptions = gh_pop_get_descendants(gs[[1]], reacttsneparent$data, path = 1)
    shiny::showModal(shiny::modalDialog(
      "Select at least 2-8 populations.",
      checkboxGroupInput("tsnepopIDs", label = "", choices = popoptions, selected = reacttsnepops$data),
      footer = list(actionButton(inputId = "tsnesettingpopIDs", label = "Ok", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "m"))
    if(length(input$tsnepopIDs) == 0){
      shinyjs::disable("tsnesettingpopIDs")
    } else{
      if(length(input$tsnepopIDs) < 2 && length(input$tsnepopIDs) > 8){
        shinyjs::disable("tsnesettingpopIDs")
      } else{
        shinyjs::enable("tsnesettingpopIDs") } } })

  observeEvent(input$tsnepopIDs,{
    if(input$tsnepopIDs[1] == ""){
      shinyjs::disable("tsnesettingpopIDs")
    } else{
      if(length(input$tsnepopIDs) < 2 && length(input$tsnepopIDs) > 8){
        shinyjs::disable("tsnesettingpopIDs")
      } else{
        shinyjs::enable("tsnesettingpopIDs") } } })

  observeEvent(input$tsnesettingpopIDs,{
    reacttsnepops$data = input$tsnepopIDs
    shiny::removeModal() })

  output$tsneplot = renderPlot({
    tsneplotactivator$data
    input$tsnehighlight
    input$tsnegrouporsampleID
    reacttsneIDs$data
    tsnegenerator <<- function(){
      par(mar = c(4, 6, 3, 26), lwd = 2)
      col = c(colorRampPalette(c("#5E4FA2", "#3288BD"))(10), colorRampPalette(c("#3288BD", "#ABDDA4"))(10), rep("#ABDDA4", 5), colorRampPalette(c("#ABDDA4", "#E6F598"))(10), colorRampPalette(c("#E6F598", "#D53E4F"))(10))
      currententiretSNE = entiretSNE
      Xlim = c(min(currententiretSNE[,1])*1.1, max(currententiretSNE[,1])*1.1)
      Ylim = c(min(currententiretSNE[,2])*1.1, max(currententiretSNE[,2])*1.1)
      if(isolate(input$tsnemode) == "Overlay Groups or Samples"){
        if(isolate(input$tsnegrouporsample) == "All"){
          shinyjs::enable("savetsneplot")
          plot(currententiretSNE[,1], currententiretSNE[,2], col = "black", xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = input$tsnedotsize/10, lwd = 0, xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
          tsneaxisticks(Xlim, Ylim)
        } else{
          if(!is.null(reacttsneIDs$data)){
            if(isolate(input$tsnegrouporsample) == "Group"){
              interval = nrow(currententiretSNE)/length(concatsamples)
              listofgrouporsamples = names(tsnelistofgroups[as.integer(reacttsneIDs$data)])
            } else if(isolate(input$tsnegrouporsample) == "Sample"){
              interval = nrow(currententiretSNE)/length(unlist(concatsamples))
              listofgrouporsamples = names(tsnelistofsamples[as.integer(reacttsneIDs$data)]) }
            dfrows = list(c(interval*as.numeric(isolate(input$tsneIDs)[1]) - interval + 1, interval*as.numeric(isolate(input$tsneIDs)[1])), c(interval*as.numeric(isolate(input$tsneIDs)[2]) - interval + 1, interval*as.numeric(isolate(input$tsneIDs)[2])))
            currenttSNE = list(currententiretSNE[dfrows[[1]][1]:dfrows[[1]][2],], currententiretSNE[dfrows[[2]][1]:dfrows[[2]][2],])
            dffinal = data.frame()
            colorlist = c("black", "red")
            for(i in seq(isolate(input$tsneIDs))){
              dffinal = rbind(dffinal, data.frame("X" = currenttSNE[[i]][,1], "Y" = currenttSNE[[i]][,2], "color" = colorlist[i])) }
            set.seed(5)
            shuffleddf = dffinal[sample(nrow(dffinal)), ]
            shinyjs::enable("savetsneplot")
            plot(shuffleddf[[1]], shuffleddf[[2]], col = shuffleddf[[3]], xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = input$tsnedotsize/10, lwd = 0, xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
            tsneaxisticks(Xlim, Ylim)
            templegend = legend("topright", bty = "n", inset = c(-0.25, -0.025), pch = 15, lty = 0, lwd = 2, cex = 1.3, pt.cex = 4, y.intersp = 1.6, text.width = strwidth("100"), pt.bg = colorlist, col = colorlist, legend = rep(NA, length(listofgrouporsamples)), xpd = TRUE)
            multiplier = 1.2
            text(strwidth(listofgrouporsamples)*multiplier + templegend$text$x*1.01, templegend$text$y, listofgrouporsamples, pos = 2, xpd = NA, font = 2, cex = multiplier) } }
        shinyjs::hide("showingtsneevents")
      } else{
        if(isolate(input$tsnegrouporsample) == "All"){
          dfrows = c(1, nrow(currententiretSNE))
          currenttSNE = currententiretSNE
        } else{
          if(isolate(input$tsnegrouporsample) == "Group"){
            interval = nrow(currententiretSNE)/length(concatsamples)
          } else if(isolate(input$tsnegrouporsample) == "Sample"){
            interval = nrow(currententiretSNE)/length(unlist(concatsamples)) }
          dfrows = c(interval*as.numeric(isolate(input$tsnegrouporsampleID)) - interval + 1, interval*as.numeric(isolate(input$tsnegrouporsampleID)))
          currenttSNE = currententiretSNE[dfrows[1]:dfrows[2],] }
        if(isolate(input$tsnemode) == "Heatmap"){
          dataf = list()
          for(i in seq(ncol(concat))){
            dataf[[i]] = concat[,i] }
          dataf = data.frame(unlist(dataf))
          dataf = cbind2(dataf, 1:nrow(dataf))
          colnames(dataf) = c("value", "ID")
          dataf = dataf %>% arrange(value)
          dataf = cbind2(dataf, colorRampPalette(col)(nrow(dataf)))
          dataf = dataf %>% arrange(ID)
          colnames(dataf)[3] = "col"
          legendrange = c(min(dataf[,1]), max(dataf[,1]))
          if(input$tsnehighlight == ""){
            highlight = availableparameters[[1]]
          } else{
            highlight = input$tsnehighlight }
          dataf = dataf[((nrow(concat)*(which(colnames(concat) == highlight) - 1)) + 1):(nrow(concat)*which(colnames(concat) == highlight)),]
          dataf["ID"] = rownames(dataf) = 1:nrow(dataf)
          shinyjs::enable("savetsneplot")
          plot(currenttSNE, xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = input$tsnedotsize/10, lwd = 0, col = dataf$col[dfrows[1]:dfrows[2]], xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
          tsneaxisticks(Xlim, Ylim)
          legendcol = colorRampPalette(col)(80)
          barXcoords = seq(Xlim[1], Xlim[2], length.out = length(legendcol))
          barposition = 1.15
          for(i in seq(legendcol)){
            points(barXcoords[i]*0.5, Ylim[2]*barposition, pch = 15, col = legendcol[i], xpd = TRUE, cex = 2) }
          text(barXcoords[1]*0.6, Ylim[2]*barposition, xpd = TRUE, labels = format(round(legendrange[1], 2), nsmall = 2), adj = 1, cex = 1.1)
          text(barXcoords[length(barXcoords)]*0.6, Ylim[2]*barposition, xpd = TRUE, labels = format(round(legendrange[2], 2), nsmall = 2), adj = 0, cex = 1.1)
          shinyjs::show("showingtsneevents")
          output$showingtsneevents = renderText({
            paste0("Plotted events: ", format((as.integer((dfrows[2] - dfrows[1]) + 1)), big.mark = ",")) })
        } else if(isolate(input$tsnemode) == "Overlay Populations" && !is.null(reacttsnepops$data)){
          groupsampindices = seq(dfrows[1], dfrows[2])
          currenttSNE = list()
          for(i in seq(reacttsnepops$data)){
            currenttSNE[[i]] = currententiretSNE[intersect(groupsampindices, which(concat[,reacttsnepops$data[i]] == 1)),] }
          dffinal = data.frame()
          colorlist = c("red2", "green4", "blue3", "orange", "black", "wheat3", "brown", "orchid3", "salmon")
          for(i in seq(reacttsnepops$data)){
            dffinal = rbind(dffinal, data.frame("X" = currenttSNE[[i]][,1], "Y" = currenttSNE[[i]][,2], "color" = colorlist[i])) }
          set.seed(6)
          shuffleddf = dffinal[sample(nrow(dffinal)), ]
          shinyjs::enable("savetsneplot")
          plot(shuffleddf[[1]], shuffleddf[[2]], col = shuffleddf[[3]], xaxt = "n", yaxt = "n", ann = FALSE, pch = 20, cex = input$tsnedotsize/10, lwd = 0, xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i")
          tsneaxisticks(Xlim, Ylim)
          templegend = legend("topright", bty = "n", inset = c(-0.25, -0.025), pch = 15, lty = 0, lwd = 2, cex = 1.3, pt.cex = 4, y.intersp = 1.6, text.width = strwidth("100"), pt.bg = colorlist[1:length(reacttsnepops$data)], col = colorlist[1:length(reacttsnepops$data)], legend = rep(NA, length(reacttsnepops$data)), xpd = TRUE)
          multiplier = 1.2
          text(strwidth(reacttsnepops$data)*multiplier + templegend$text$x*1.01, templegend$text$y, reacttsnepops$data, pos = 2, xpd = NA, font = 2, cex = multiplier)
          shinyjs::hide("showingtsneevents") } } }
    if(!is.null(entiretSNE)){
      shinyjs::disable("savetsneplot")
      tsnegenerator()
      shinyjs::enable("tsnegenerate")
    } else{
      shinyjs::disable("savetsneplot") }
  }, height = 440, width = 800)

  output$savetsneplot = downloadHandler(
    filename = "t-SNE plot.png",
    content = function(file){
      png(file, units = "in", height = 6, width = 11.02, res = 300)
      tsnegenerator()
      dev.off() } )

  #Results----
  ##left----
  observeEvent(input$rsstat, {
    if(isolate(input$rsstat) == "Freq. of parent" || isolate(input$rsstat) == "Freq. of total" || isolate(input$rsstat) == "Count"){
      shinyjs::disable("rsparent")
      shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
      shinyjs::disable("rsparameter")
    } else{
      if(isolate(input$rsstat) == "Freq. of..."){
        shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
        shinyjs::disable("rsparameter")
      } else{
        shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
        shinyjs::enable("rsparameter") } } })

  observeEvent(c(input$rsstat, input$rsparent, input$rspop, input$rsparameter), {
    if(isolate(input$rspop) != ""){
      if(isolate(input$rsstat) == "Freq. of parent" || isolate(input$rsstat) == "Freq. of total" || isolate(input$rsstat) == "Count"){
        shinyjs::enable("addresult")
      } else{
        if(isolate(input$rsparent) != "" || isolate(input$rsparameter) != ""){
          shinyjs::enable("addresult")
        } else{
          shinyjs::disable("addresult") } }
    } else{
      shinyjs::disable("addresult") } })

  observeEvent(c(input$rsstat, input$rspop), {
    if(!is.null(gs)){
      fullpath = flowWorkspace::gs_get_pop_paths(gs)[which(flowWorkspace::gs_get_pop_paths(gs, path = 1) == input$rspop)]
      possibleparents = strsplit(fullpath, split = "/")
      if(length(possibleparents) > 0){
        possibleparents = possibleparents[[1]][1:length(possibleparents[[1]]) - 1] }
      if(isolate(input$rsstat) == "Freq. of..." && isolate(input$rspop) != ""){
        shiny::updateSelectInput(inputId = "rsparent", choices = c(possibleparents))
        shinyjs::enable("rsparent")
      } else{
        shiny::updateSelectInput(inputId = "rsparent", selected = "")
        shinyjs::disable("rsparent") } } })

  observeEvent(input$addresult, {
    shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
    if(is.na(results[[1]][1])){
      resultID = 1
    } else{
      resultID = ncol(results) + 1 }
    popstats = flowWorkspace::gs_pop_get_count_fast(gs, subpopulations = input$rspop)
    if(input$rsstat == "Freq. of parent"){
      results[[resultID]] <<- sprintf("%.1f", popstats[[4]]*100/popstats[[5]])[onlysampleIDs]
    } else if(input$rsstat == "Freq. of total"){
      results[[resultID]] <<- sprintf("%.1f", flowWorkspace::gs_pop_get_count_fast(gs, "freq", subpopulations = input$rspop)[,4]$Frequency*100)[onlysampleIDs]
    } else if(input$rsstat == "Count"){
      results[[resultID]] <<- popstats[[4]][onlysampleIDs]
    } else if(input$rsstat == "Freq. of..."){
      parentstats = flowWorkspace::gs_pop_get_count_fast(gs, subpopulations = input$rsparent)
      results[[resultID]] <<- sprintf("%.1f", popstats[[4]]*100/parentstats[[4]])[onlysampleIDs]
    } else if(input$rsstat == "Median"){
      results[[resultID]] <<- sprintf("%.1f", flowWorkspace::gs_pop_get_stats(gs, input$rspop, pop.MFI)[[which(input$rsparameter == fluochannelsfull) + 2]]) }
    if(input$rsstat == "Median"){
      name1 = "MFI"
    } else{
      name1 = input$rsstat }
    if(input$rsparent != ""){
      name2 = paste0(input$rsparent)
    } else{
      name2 = "" }
    if(input$rsparameter != ""){
      name4 = paste0(" (", input$rsparameter, ")")
    } else{
      name4 = "" }
    colnames(results)[resultID] <<- paste0(name1, name2, "/", input$rspop, name4) })

  ##main----
  output$results = renderTable(rownames = TRUE, spacing = "xs", {
    input$addresult
    input$deletemodal
    input$tabs
    if(!is.na(results[,1][[1]])){
      shinyjs::enable("exporttable")
      shinyjs::enable("editresult") }
    allpops = flowWorkspace::gs_get_pop_paths(gs, path = 1)
    allpops = allpops[allpops != allpops[1]]
    shiny::updateSelectInput(inputId = "rspop", choices = c("", allpops))
    results})

  observeEvent(input$editresult, {
    shiny::showModal(shiny::modalDialog(
      "Select the rows to be deleted.",
      checkboxGroupInput("checkdelete", label = "", choices = setNames(as.list(seq(ncol(results))), colnames(results)) ),
      footer = list(actionButton(inputId = "deletemodal", label = "Delete", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color:white; border-color:black")),
      easyClose = TRUE,
      size = "s")) })

  observeEvent(input$deletemodal, {
    results <<- results[-as.integer(input$checkdelete)]
    shiny::removeModal()
    if(ncol(results) == 0){
      results[1] <<- NA
      colnames(results) <<- NA
      shinyjs::disable("exporttable")
      shinyjs::disable("editresult") } })

  output$exporttable = downloadHandler(
    filename = "Results.csv",
    content = function(file){
      write.csv(results, file) } )

  #Between Tabs----
  observeEvent(input$tabs,{
    if(input$tabs == "comptab"){
      uncompgsfortable <<- flowWorkspace::gs_clone(uncompgs)
      for(i in seq(fluochannels)){
        uncompgsfortable <<- transform(uncompgsfortable, transformerList(fluochannels[i], flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3]))) }
      shinyjs::hide("savematrix")
      shinyjs::hide("cancelmatrix")
      shinyjs::enable("creatematrix")
      reactreadnonly$data = TRUE
      hTable$data = as.data.frame(compdfs["Cytometer-defined"][[1]]*100)
      output$comp = rhandsontable::renderRHandsontable({
        #BUG----
        #the %>% comes from rhandsontable::
        visualtable = rhandsontable::rhandsontable(compdfs[input$previewmatrix][[1]]*100, rowHeaderWidth = 100, stretchH = "all", readOnly = reactreadnonly$data) %>%
          rhandsontable::hot_validate_numeric(cols = seq(fluochannels)) %>%
          rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE) %>%
          rhandsontable::hot_cols(format = "0", renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                 Handsontable.renderers.NumericRenderer.apply(this, arguments);
                 heatmap = ['#FFFFCC66', '#FFEDA066', '#FED97666', '#FEB24C66', '#FD8D3C66', '#FC4E2A66', '#E31A1C66', '#BD002666', '#80002666']
                 if (row == col) {
                   td.style.background = 'lightgrey';
                 } else if (Math.abs(value) > 88) {
                   td.style.background = heatmap[8];
                 } else if (Math.abs(value) > 77) {
                   td.style.background = heatmap[7];
                 } else if (Math.abs(value) > 66) {
                   td.style.background = heatmap[6];
                 } else if (Math.abs(value) > 55) {
                   td.style.background = heatmap[5];
                 } else if (Math.abs(value) > 44) {
                   td.style.background = heatmap[4];
                 } else if (Math.abs(value) > 33) {
                   td.style.background = heatmap[3];
                 } else if (Math.abs(value) > 22) {
                   td.style.background = heatmap[2];
                 } else if (Math.abs(value) > 11) {
                   td.style.background = heatmap[1];
                 } else if (Math.abs(value) > 0) {
                   td.style.background = heatmap[0];
                 } else if (Math.abs(value) <= 0) {
                   td.style.background = 'white'} }")
        for(i in seq(fluochannels)){
          cell = list(row = i - 1, col = i - 1, readOnly = TRUE)
          visualtable$x$cell = c(visualtable$x$cell, list(cell)) }
        htmlwidgets::onRender(visualtable, change_hook) })
    } else if(input$tabs == "proliftab"){
      poppaths = flowWorkspace::gs_get_pop_paths(gs, path = 1)
      if("undivided" %in% poppaths){
        shinyjs::hide("prolifstart")
        shinyjs::show("prolifplot")
      } else{
        shinyjs::hide("prolifplot")
        shinyjs::show("prolifstart") } } })

  #File----
  shinyFiles::shinyFileChoose(input, "BCytofileload", roots = c(Home = fs::path_home()), session = session, restrictions = system.file(package = "base"), filetypes = "RData")
  shinyFiles::shinyDirChoose(input, "directory", roots = c(Home = fs::path_home()), session = session, restrictions = system.file(package = "base"), allowDirCreate = FALSE, filetypes = "fcs")

  observeEvent(input$directory,{
    if(!is.integer(input$directory[1])){
      reactdir$data = input$directory } })

  observeEvent(reactdir$data,{
    shiny::withProgress(message = "Please wait...", detail = "", value = 0, max = 100, {
      tempfilepath <<- shinyFiles::parseDirPath(c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()), reactdir$data)
      tempfilelist <<- list.files(tempfilepath, pattern = ".fcs")
      if(length(tempfilelist) == 0){
        shinyjs::alert("No FCS file was detected in the selected folder.")
      } else{
        setwd(tempfilepath)
        validityout = list()
        for(i in seq(tempfilelist)){
          validityout[[i]] = tryCatch({flowCore::keyword(flowWorkspace::load_cytoframe_from_fcs(tempfilelist[i]))
          }, trycatchfail = function(x){ } ) }
        if(all(sapply(validityout, function(x) length(x) > 0))){
          channellengths = sapply(tempfilelist, function(x) length(colnames(flowWorkspace::load_cytoframe_from_fcs(x))))
          if(length(unique(channellengths)) == 1){
            listoftubenames = list()
            listoftubes = list()
            for(i in seq(tempfilelist)){
              listoftubes[[i]] = validityout[[i]]$ORIGINALITY
              listoftubenames[[i]] = validityout[[i]]$`TUBE NAME` }
            if(length(listoftubenames) > 0){
              tempcompcontrolIDs <<- list()
              for(i in seq(tempfilelist)){
                tubename = strsplit(listoftubenames[[i]], " ")[[1]]
                if(length(grep("Control", tubename)) > 0){
                  if(length(grep("Stained", tubename)) > 0 || length(grep("Unstained", tubename)) > 0){
                    tempcompcontrolIDs[[length(tempcompcontrolIDs) + 1]] = i } } }
              tempcompcontrolIDs = unlist(tempcompcontrolIDs)
              tempcompcontrolIDs <<- tempcompcontrolIDs
              finalselection <<- list(as.list(seq(tempfilelist)), as.list(seq(tempfilelist)))
              if(!is.null(tempcompcontrolIDs)){
                temponlysampleIDs <<- seq(tempfilelist)[-tempcompcontrolIDs]
                finalselection[[1]][tempcompcontrolIDs] <<- TRUE
                finalselection[[1]][-tempcompcontrolIDs] <<- FALSE
                finalselection[[2]][temponlysampleIDs] <<- TRUE
                finalselection[[2]][-temponlysampleIDs] <<- FALSE
                if(length(temponlysampleIDs) == 0){
                  finalselection[[2]][tempcompcontrolIDs] = FALSE }
              } else{
                finalselection[[1]] <<- lapply(finalselection[[1]], function(x) FALSE)
                finalselection[[2]] <<- lapply(finalselection[[2]], function(x) TRUE) }
              samplecolum = list()
              for(i in seq(tempfilelist)){
                samplecolum[[i]] = list(br(), br(), tempfilelist[i]) }
              groupcolums = list(list(), list())
              for(j in seq(2)){
                for(i in seq(tempfilelist)){
                  groupcolums[[j]][[i]] = checkboxInput(inputId = paste0("column", j, "sample", i), label = "", value = finalselection[[j]][[i]]) } }
              checkboxcolumns = list(column(width = 4, align = "center", style = "margin-top: 10px;", list(strong("Compensation Controls"), groupcolums[[1]])))
              checkboxcolumns[[2]] = column(width = 2, align = "center", style = "margin-top: 10px;", list(strong("Samples"), groupcolums[[2]]))
              shiny::showModal(shiny::modalDialog(
                strong("Please indicate/confirm compensation controls and sample files."),
                br(),
                fluidRow(column(width = 6, align = "right", list(samplecolum)),
                         checkboxcolumns ),
                footer = list(actionButton(inputId = "okselectcompandsamp", label = "Ok", style="color: white; background-color:black; border-color:black"),actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color: white; border-color: black")),
                easyClose = FALSE,
                size = "l"))
            } else{
              if(length(listoftubes) > 0){
                shinyjs::alert("Loading failed... At least one FCS file had its integrity altered. BCyto supports only original unmodified FCS files.")
              } else{
                shinyjs::alert("Loading failed... At least one FCS file is of an unsupported version or is corrupted.") } }
          } else{
            shinyjs::alert("Loading failed... This software requires all FCS files to be consistent with each other in the number of channels.") }
        } else{
          shinyjs::alert("Loading failed... At least one FCS file seems to be corrupted.") } }
      reactdir$data = NULL }) })

  observeEvent(input$okselectcompandsamp,{
    for(j in seq(2)){
      for(i in seq(tempfilelist)){
        finalselection[[j]][[i]] <<- input[[paste0("column", j, "sample", i)]] } }
    displayalert = NULL
    for(i in seq(tempfilelist)){
      if(finalselection[[1]][[i]] == FALSE && finalselection[[2]][[i]] == FALSE){
        displayalert = "At least one file has not been checked."
      } else if(finalselection[[1]][[i]] == TRUE && finalselection[[2]][[i]] == TRUE){
        displayalert = "At least one file has been checked twice." } }
    if(all(sapply(finalselection[[2]], function(x) x == FALSE))){
      displayalert = "Only compensation control files were loaded. This software requires at least a sample file to be loaded." }
    if(is.null(displayalert)){
      shiny::withProgress(message = "Please wait...", detail = "", value = 0, max = 100, {
        shinyjs::hide("okselectcompandsamp")
        shinyjs::hide("cancelmodal")
        filepath <<- tempfilepath
        filelist <<- tempfilelist
        compcontrolIDs <<- seq(filelist)[unlist(finalselection[[1]])]
        onlysampleIDs <<- seq(filelist)[unlist(finalselection[[2]])]
        initializing(filelist)
        plotactivator$data = plotactivator$data + 1
        hierarchyactivator$data = hierarchyactivator$data + 1
        choices = names(compdfs)
        names(choices) = choices
        allnames = c()
        for(i in seq(choices)){
          if(choices[i] == appliedmatrix){
            allnames[i] = paste(appliedmatrix, "(applied)")
          } else{
            allnames[i] = choices[i] } }
        names(choices) = allnames
        shiny::updateSelectInput(inputId = "previewmatrix", choices = choices, selected = appliedmatrix)
        shiny::updateSelectInput(inputId = "bgtype", selected = "")
        shiny::updateSelectInput(inputId = "ovtype", selected = "")
        shinyjs::js$disableTab("ancestrytab")
        shinyjs::js$disableTab("overlaytab")
        shinyjs::js$disableTab("proliftab")
        shinyjs::js$disableTab("tsnetab")
        shinyjs::js$disableTab("resulttab")
        shinyjs::disable("editresult")
        shinyjs::disable("exporttable")
        shinyjs::hide("tsnemode")
        shinyjs::hide("tsnedotsize")
        shinyjs::hide("tsnehighlight")
        shinyjs::hide("tsnegrouporsample")
        shinyjs::hide("tsnegrouporsampleID")
        shinyjs::hide("showingtsneevents")
        shinyjs::hide("savetsneplot")
        shinyjs::enable("savefile")
        shinyjs::enable("tsnegenerate")
        shinyjs::show("tsnegroups")
        shinyjs::show("tsnesamples")
        shinyjs::show("tsneparent")
        shinyjs::show("tsneparameters")
        shinyjs::show("tsneevents")
        betweenpeaks <<- refpeaks <<- entiretSNE <<- tsnelistofgroups <<- tsnelistofsamples <<- concat <<- concatsamples <<- samples <<- NULL
        reacttsneparent$data = NULL
        reacttsnepar$data = NULL
        shiny::updateSelectInput(inputId = "tsneevents", selected = 0)
        reactgroupnames$data = NULL
        reacttsneIDs$data = NULL
        reacttsnepops$data = NULL
        reactovsamples$data = NULL
        tsneplotactivator$data = tsneplotactivator$data + 1
        if(length(fluochannels) >= 10){
          reactheight$data = 850
        } else if(length(fluochannels) < 10 && length(fluochannels) > 4){
          reactheight$data = length(fluochannels)*85
        } else if(length(fluochannels) <= 4){
          reactheight$data = 350 }
        shinyjs::alert("Files have been loaded successfully.")
        shiny::removeModal() })
    } else{
      setwd(filepath)
      shinyjs::alert(displayalert) } })

  observeEvent(input$loadtestdata,{
    shiny::showModal(shiny::modalDialog(
      strong("Select the data set to be downloaded."),
      br(),
      radioButtons("radiotestdata", label = "", choiceValues = c("test-data-1", "test-data-2"), choiceNames = c("Cell maturation assay - 24 mb (FlowRepository ID FR-FCM-Z4LZ)", "Proliferation assay - 18.5 mb (FlowRepository ID FR-FCM-Z4LY)"), width = "100%"),
      footer = list(actionButton(inputId = "downloadtestdata", label = "Download", style="color: white; background-color:black; border-color:black"), actionButton(inputId = "cancelmodal", label = "Cancel", style="color: black; background-color: white; border-color: black")),
      easyClose = FALSE,
      size = "m")) })

  observeEvent(input$downloadtestdata,{
    shinyjs::hide("downloadtestdata")
    shinyjs::hide("cancelmodal")
    if(input$radiotestdata == "test-data-1"){
      filelist <<- c("Compensation Controls_APC Stained Control_013.fcs", "Compensation Controls_FITC Stained Control_011.fcs", "Compensation Controls_PE Stained Control_014.fcs", "Compensation Controls_PE-Cy7 Stained Control_015.fcs", "Compensation Controls_Unstained Control_010.fcs", "Compensation Controls_V 450 Stained Control_012.fcs", "Specimen_001_FMO-CD11c_003.fcs", "Specimen_001_FMO-CD86_005.fcs", "Specimen_001_FMO-Gr1_002.fcs", "Specimen_001_FMO-MHCII_004.fcs", "Specimen_001_FMO-Viability_001.fcs", "Specimen_001_iDC_001_006.fcs", "Specimen_001_iDC_002_007.fcs", "Specimen_001_mDC_001_008.fcs", "Specimen_001_mDC_002_009.fcs")
      compcontrolIDs <<- c(1:6)
      onlysampleIDs <<- c(7:15)
      datasetfile = "dataset-1.RData"
    } else{
      filelist <<- c("Compensation Controls_APC-Cy7 Stained Control_080.fcs", "Compensation Controls_FITC Stained Control_078.fcs", "Compensation Controls_PE Stained Control_081.fcs", "Compensation Controls_Unstained Control_077.fcs", "Compensation Controls_V 450 Stained Control_079.fcs", "Specimen_001_0-00ug-mL_002_030.fcs", "Specimen_001_0-01ug-mL_003_034.fcs", "Specimen_001_0-05ug-mL_002_036.fcs", "Specimen_001_1-00ug-mL_001_044.fcs", "Specimen_001_5-00ug-mL_002_048.fcs", "Specimen_001_FMO-CD4_003.fcs", "Specimen_001_FMO-CD44_004.fcs", "Specimen_001_FMO-CFSE_001.fcs", "Specimen_001_FMO-Viability_002.fcs")
      compcontrolIDs <<- c(1:5)
      onlysampleIDs <<- c(6:14)
      datasetfile = "dataset-2.RData" }
    filepath <<- tempdir()
    if(all(sapply(filelist, function(x) length(grep(x, list.files(filepath))) > 0)) == FALSE){
      adjustedfilelist = list()
      for(i in seq(filelist)){
        if(grepl(" ", filelist[i])){
          adjustedfilelist[[i]] = gsub(" ", "%20", filelist[i])
        } else{
          adjustedfilelist[[i]] = filelist[i] } }
      adjustedfilelist = unlist(adjustedfilelist)
      dat = data.frame(x = numeric(0), y = numeric(0))
      shiny::withProgress(message = "Downloading", detail = "file 0", value = 0, max = (length(filelist) + 1), {
        for(i in seq(filelist)){
          dat = rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
          shiny::incProgress(1, detail = paste("file", i))
          Sys.sleep(0.1)
          download.file(paste0("https://github.com/BonilhaCaio/", input$radiotestdata, "/raw/main/", adjustedfilelist[i]), destfile = file.path(filepath, filelist[i])) }
        download.file(paste0("https://github.com/BonilhaCaio/", input$radiotestdata, "/raw/main/", datasetfile), destfile = file.path(filepath, datasetfile)) }) }
    reactBCytofile$data = datasetfile
    if(input$radiotestdata == "test-data-1"){
      shinyjs::js$disableTab("proliftab")
    } else{
      shinyjs::js$enableTab("proliftab") }
    shinyjs::js$enableTab("ancestrytab")
    shinyjs::js$enableTab("overlaytab")
    shinyjs::js$enableTab("tsnetab")
    shinyjs::js$enableTab("resulttab") })

  observeEvent(input$BCytofileload,{
    if(!is.integer(input$BCytofileload[1])){
      reactBCytofile$data = input$BCytofileload } })

  observeEvent(reactBCytofile$data,{
    if(reactBCytofile$data[[1]] == "dataset-1.RData" || reactBCytofile$data[[1]] == "dataset-2.RData"){
      booleantest = TRUE
    } else{
      tempfilepath = shinyFiles::parseFilePaths(c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()), reactBCytofile$data)
      tempstrsplit = strsplit(tempfilepath[[4]][[1]], "")[[1]]
      filepath <<- paste(tempstrsplit[-(length(tempstrsplit) - length(strsplit(tempfilepath[[1]], "")[[1]])):-length(tempstrsplit)], collapse = "")
      load(paste0(filepath, "/", tempfilepath[[1]]), testenvironment <- new.env())
      booleantest = all(sapply(testenvironment$filelist, function(x) length(grep(x, list.files(filepath))) > 0)) }
    if(booleantest){
      loadingfile <<- TRUE
      shiny::withProgress(message = "Please wait...", detail = "", value = 0, max = 100, {
        setwd(filepath)
        if(reactBCytofile$data[[1]] == "dataset-1.RData" || reactBCytofile$data[[1]] == "dataset-2.RData"){
          load(paste0(filepath, "/", reactBCytofile$data))
        } else{
          load(paste0(filepath, "/", tempfilepath[[1]])) }
        uncompgs <<- flowWorkspace::GatingSet(flowWorkspace::load_cytoset_from_fcs(filelist))
        gs <<- flowWorkspace::gs_clone(uncompgs)
        flowCore::compensate(gs, compdfs[appliedmatrix][[1]])
        allchannels <<- colnames(gs)

        print(allchannels)

        fluochannels <<- allchannels[grepl("SC", allchannels) != TRUE][allchannels[grepl("SC", allchannels) != TRUE] != "Time"]
        hTable <<- shiny::reactiveValues(data = as.data.frame(compdfs[[1]]*100))
        for(i in seq(fluochannels)){
          gs = transform(gs, transformerList(fluochannels[i], flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=currentcustomaxis[[i]][1], neg=currentcustomaxis[[i]][2], widthBasis=currentcustomaxis[[i]][3]))) }
        for(i in seq(popgates)){
          flowWorkspace::gs_pop_add(gs, popgates[[i]], parent = popparents[i]) }
        flowWorkspace::recompute(gs)
        reactshowaxis$data = plotshowaxis
        reactaxisfont$data = plotaxisfont
        reactshowgatename$data = plotshowgatename
        reactgatefont$data = plotgatefont
        choices = names(compdfs)
        names(choices) = choices
        allnames = c()
        for(i in seq(choices)){
          if(choices[i] == appliedmatrix){
            allnames[i] = paste(appliedmatrix, "(applied)")
          } else{
            allnames[i] = choices[i] } }
        names(choices) = allnames
        shiny::updateSelectInput(inputId = "previewmatrix", choices = choices, selected = appliedmatrix)
        pops = flowWorkspace::gs_get_pop_paths(gs, path = 1)
        shiny::updateSelectInput(inputId = "bgsample", choices = c("", filelist[onlysampleIDs]), selected = bgsample)
        shiny::updateSelectInput(inputId = "bgtype", selected = bgtype)
        shiny::updateSelectInput(inputId = "bgpop", choices = c("", pops[-1]), selected = bgpop)
        showovplot$data = ovboolean
        shiny::updateSelectInput(inputId = "ovtype", selected = ovtype)
        shiny::updateSelectInput(inputId = "ovtone", selected = ovtone)
        shiny::updateSelectInput(inputId = "ovchannelY", choices = c(ch), selected = ovchannelY)
        shiny::updateSelectInput(inputId = "ovchannelX", choices = c(ch), selected = ovchannelX)
        shiny::updateSelectInput(inputId = "ovparent", choices = pops, selected = ovparent)
        reactovsamples$data = ovsamples
        reactovshowaxis$data = ovshowaxis
        reactovaxisfont$data = ovaxisfont
        prolifready$data = prolifboolean
        shiny::updateSelectInput(inputId = "prolifsamplebutt", choices = filelist[onlysampleIDs], selected = prolifsamplebutt)
        if(prolifboolean == TRUE){
          updateCheckboxInput(inputId = "proliflabel", value = proliflabel)
          updateCheckboxInput(inputId = "prolifgrid", value = prolifgrid) }
        reacttsneparent$data = tsneparent
        reacttsnepar$data = tsneparameters
        if(is.null(tsneevents)){
          tsneevents = 0 }
        shiny::updateSelectInput(inputId = "tsneevents", selected = tsneevents)
        reactgroupnames$data = tsnegroupnames
        if(length(tsnelistofgroups) < 2 && length(tsnelistofsamples) < 2){
          shiny::updateSelectInput(inputId = "tsnemode", choices = c("Heatmap", "Overlay Populations"), selected = tsnemode)
        } else{
          shiny::updateSelectInput(inputId = "tsnemode", choices = c("Heatmap", "Overlay Groups or Samples", "Overlay Populations"), selected = tsnemode) }
        shiny::updateSelectInput(inputId = "tsnehighlight", selected = tsnehighlight)
        shiny::updateSelectInput(inputId = "tsnegrouporsample", selected = tsnegrouporsample)
        shiny::updateSelectInput(inputId = "tsnegrouporsampleID", selected = tsnegrouporsampleID)
        reacttsneIDs$data = tsnegrouporsampleIDs
        reacttsnepops$data = tsnepops
        updateSliderInput(inputId = "tsnedotsize", value = tsnedotsize)
        if(!is.null(entiretSNE)){
          shinyjs::hide("tsnegroups")
          shinyjs::hide("tsnesamples")
          shinyjs::hide("tsneparent")
          shinyjs::hide("tsneparameters")
          shinyjs::hide("tsneevents")
          shinyjs::show("tsnemode")
          shinyjs::show("tsnedotsize")
          shinyjs::show("tsnehighlight")
          shinyjs::show("tsnegrouporsample")
          shinyjs::show("tsnegrouporsampleID")
          shinyjs::show("showingtsneevents")
          shinyjs::show("savetsneplot")
          shinyjs::enable("tsnegenerate")
          tempchannels = unlist(ch)[unlist(ch) != "Time"]
          availableparameters <<- list()
          for(i in seq(tsneparameters)){
            availableparameters[[i]] <<- which(tempchannels == tsneparameters[i]) }
          availableparameters <<- tempchannels[unlist(availableparameters)]
          shiny::updateSelectInput(inputId = "tsnehighlight", choices = c(availableparameters))
          tsneplotactivator$data = tsneplotactivator$data + 1
        } else{
          shinyjs::hide("tsnemode")
          shinyjs::hide("tsnedotsize")
          shinyjs::hide("tsnehighlight")
          shinyjs::hide("tsnegrouporsample")
          shinyjs::hide("tsnegrouporsampleID")
          shinyjs::hide("showingtsneevents")
          shinyjs::hide("savetsneplot")
          shinyjs::show("tsnegroups")
          shinyjs::show("tsnesamples")
          shinyjs::show("tsneparent")
          shinyjs::show("tsneparameters")
          shinyjs::show("tsneevents")
          shinyjs::enable("tsnegenerate") } })
      filelist <<- filelist
      onlysampleIDs <<- onlysampleIDs
      currentcustomaxis <<- currentcustomaxis
      compdfs <<- compdfs
      appliedmatrix <<- appliedmatrix
      betweenpeaks <<- betweenpeaks
      refpeaks <<- refpeaks
      adjustedlabels <<- adjustedlabels
      prolifchannel <<- prolifchannel
      prolifparent <<- prolifparent
      entiretSNE <<- entiretSNE
      ch <<- ch
      tsnelistofgroups <<- tsnelistofgroups
      tsnelistofsamples <<- tsnelistofsamples
      concat <<- concat
      concatsamples <<- concatsamples
      results <<- results
      plotactivator$data = plotactivator$data + 1
      X <<- allchannels[grep("FS", allchannels)[1]]
      Y <<- allchannels[grep("SS", allchannels)[1]]
      fluochannelsfull <<- names(ch[grepl("SC", ch) != TRUE][ch[grepl("SC", ch) != TRUE] != "Time"])
      shiny::updateSelectInput(inputId = "sample", choices = c("", filelist[onlysampleIDs]), selected = filelist[onlysampleIDs][1])
      shiny::updateSelectInput(inputId = "type", selected = "Pseudocolor")
      shiny::updateSelectInput(inputId = "previewsample", choices = c("", filelist), selected = filelist[onlysampleIDs][1])
      shiny::updateSelectInput(inputId = "channelY", choices = c(ch), selected = Y)
      shiny::updateSelectInput(inputId = "channelX", choices = c(ch), selected = X)
      shiny::updateSelectInput(inputId = "rsparameter", choices = c("", fluochannelsfull))
      shinyjs::js$enableTab("plottab")
      shinyjs::js$enableTab("comptab")
      shinyjs::js$disableTab("ancestrytab")
      shinyjs::js$disableTab("overlaytab")
      shinyjs::js$disableTab("proliftab")
      shinyjs::js$disableTab("tsnetab")
      shinyjs::js$disableTab("resulttab")
      shiny::removeModal()
      shinyjs::alert("File has been loaded successfully.")
      shinyjs::enable("savefile")
      reactparent$data = "root"
      hierarchyactivator$data = hierarchyactivator$data + 1
      reactBCytofile$data = NULL
      if(length(fluochannels) >= 10){
        reactheight$data = 850
      } else if(length(fluochannels) < 10 && length(fluochannels) > 4){
        reactheight$data = length(fluochannels)*85
      } else if(length(fluochannels) <= 4){
        reactheight$data = 350 }
    } else{
      shinyjs::alert("Loading failed... Please certify this save file and all FCS files connected to it are in the folder.") } })

  observeEvent(input$savefile,{
    setwd(filepath)
    poppaths = flowWorkspace::gs_get_pop_paths(gs)
    popparents = flowWorkspace::gs_pop_get_count_fast(gs[1], "freq")[,3][[1]]
    popgates = list()
    for(i in seq(poppaths)[-1]){
      popgates[[i]] = flowWorkspace::gs_pop_get_gate(gs, poppaths[i])[1][[1]] }
    popgates = popgates[-1]
    plotshowaxis = reactshowaxis$data
    plotaxisfont = reactaxisfont$data
    plotshowgatename = reactshowgatename$data
    plotgatefont = reactgatefont$data
    bgsample = input$bgsample
    bgtype = input$bgtype
    bgpop = input$bgpop
    ovboolean = showovplot$data
    ovtype = input$ovtype
    ovtone = input$ovtone
    ovchannelY = input$ovchannelY
    ovchannelX = input$ovchannelX
    ovparent = input$ovparent
    ovsamples = reactovsamples$data
    ovshowaxis = reactovshowaxis$data
    ovaxisfont = reactovaxisfont$data
    prolifboolean = prolifready$data
    prolifsamplebutt = input$prolifsamplebutt
    proliflabel = input$proliflabel
    prolifgrid = input$prolifgrid
    tsneparent = reacttsneparent$data
    tsneparameters = reacttsnepar$data
    tsneevents = input$tsneevents
    tsnegroupnames = reactgroupnames$data
    tsnemode = input$tsnemode
    tsnehighlight = input$tsnehighlight
    tsnegrouporsample = input$tsnegrouporsample
    tsnegrouporsampleID = input$tsnegrouporsampleID
    tsnegrouporsampleIDs = reacttsneIDs$data
    tsnepops = reacttsnepops$data
    tsnedotsize = input$tsnedotsize
    save(filelist, onlysampleIDs, currentcustomaxis, popparents, popgates, plotshowaxis, plotaxisfont, plotshowgatename, plotgatefont, ch, #Plot
         compdfs, appliedmatrix, #Compensation
         bgsample, bgtype, bgpop, #Ancestry
         ovboolean, ovtype, ovtone, ovchannelY, ovchannelX, ovparent, ovsamples, ovshowaxis, ovaxisfont, #Overlays
         prolifchannel, prolifparent, prolifboolean, betweenpeaks, refpeaks, adjustedlabels, prolifsamplebutt, proliflabel, prolifgrid, #Proliferation
         concat, entiretSNE, concatsamples, tsneparent, tsneparameters, tsneevents, tsnegroupnames, tsnemode, tsnehighlight, tsnegrouporsample, tsnegrouporsampleID, tsnegrouporsampleIDs, tsnepops, tsnedotsize, tsnelistofgroups, tsnelistofsamples, #t-SNE
         results, #Results
         file = paste0(getwd(), "/BCytoSave.RData"))
    shinyjs::alert(paste("A BCyto file has been saved in the same folder of the FCS files:", filepath)) })

  observeEvent(input$filehelp, {
    test <<- input$directory
    shiny::showModal(shiny::modalDialog(
      title = "'File' tab",
      "To start, load your FCS files, a BCyto file or download test data.",
      br(), br(), br(),
      strong("BCyto files need to be kept inside the folder of corresponding FCS files."),
      footer = NULL,
      size = "m",
      easyClose = TRUE)) })





}








  shiny::shinyApp(ui, server) }
runBCyto()




#NOTES FOR QUESTION MARK ICON:
#For estimating the count on the Y axis on histograms, the X axis was divided in 20 bins (hist(dfX, breaks...))
#Normalization on the overlay plots (Y axis) is made by dividing the density of the points made for drawing the lines by the number of events, to allow normalization between samples with differentn number of gated events
#Lines in contour plots utilize a sampling method to select up to 5,000 events. Sampling is done by utilizing sample.int() to generate random row numbers, representing events in FACS data, from the gatingset matrix. Reproducibility is ensured by set.seed(). It is important to note that this sampling does not affect statistics generated by gates and it is only for aesthetics purpose, as the full matrix is used for that.
#In the Ancestry tab plots, the overlay is univariate, but in the Overlay tab plots, dot plots are bivariate (the two-color points are individually added in mixed orders) (https://docs.flowjo.com/flowjo/graphical-reports/graph-options-and-annotation/le-overlays/)
#Following transformation with flowWorkspace::flowWorkspace::flowjo_biexp_trans(), the ggcyto R package was used to acquired the numbers for axis ticks of fluorescent channels. A plot was built with autoplot() with the transformed gating set, and breaks were accessed from the $axis_inverse_trans slot.
#for acquiring identical plots in Flowjo, for comparison with this software, axis scales from fluorescence channels should be set to Extra Neg. Decades 0.8, Width Basis -10 and Positive Decades 4.5

#AXIS LIMIT CUSTOMIZATION equivalent/legend to put in the info (all this and channelRange=4096, maxValue=262144, widthBasis=-10):
#min 10^-3 (Negative Decades = 0), max 10^5 (Positive Decades = 4.42)
#min 10^-2 (Negative Decades = 1), max 10^5 (Positive Decades = 4.42)
#min 10^-3 (Negative Decades = 0), max 10^4 (Positive Decades = 3.42)
#min 10^-2 (Negative Decades = 1), max 10^4 (Positive Decades = 3.42)
