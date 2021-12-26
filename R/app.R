#BCyto----
#' Analyse and visualize Flow Cytometry data
#' @description The BCyto R package is an open-source project that
#' provides an user-friendly, high-performance UI for Flow Cytometry analysis.
#' @docType package
#' @name BCyto
NULL
#' Runs the BCyto app
#' @description runBCyto() is a wrapper function of
#' shinyApp(ui, server, onStart).
#' @details Use runBCyto() at the R console for initiating the BCyto Shiny app.
#'
#' @export

runBCyto <- function() {
    shinyApp(ui, server, onStart = function() {
        useShinyjs()
        extendShinyjs(text=jsCode, functions=c('disableTab','enableTab'))
        inlineCSS(css)
    })
}

#' @import autospill
#' @import data.table
#' @import graphics
#' @import shiny
#' @importClassesFrom flowCore flowFrame
#' @importClassesFrom flowCore flowSet
#' @importClassesFrom flowWorkspace GatingHierarchy
#' @importClassesFrom flowWorkspace GatingSet
#' @importClassesFrom flowWorkspace booleanFilter
#' @importClassesFrom flowWorkspace cytoframe
#' @importClassesFrom flowWorkspace cytoset
#' @importFrom MASS kde2d
#' @importFrom Rgraphviz agopen
#' @importFrom Rtsne Rtsne
#' @importFrom base64enc dataURI
#' @importFrom dplyr arrange
#' @importFrom flowCore compensate
#' @importFrom flowCore flowSet_to_list
#' @importFrom flowCore keyword
#' @importFrom flowCore polygonGate
#' @importFrom flowCore quadGate
#' @importFrom flowCore rectangleGate
#' @importFrom flowCore spillover
#' @importFrom flowWorkspace GatingSet
#' @importFrom flowWorkspace cytoframe_to_flowFrame
#' @importFrom flowWorkspace cytoset
#' @importFrom flowWorkspace cytoset_to_flowSet
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom flowWorkspace flowjo_biexp_trans
#' @importFrom flowWorkspace gh_pop_get_descendants
#' @importFrom flowWorkspace gh_pop_get_indices
#' @importFrom flowWorkspace gs_clone
#' @importFrom flowWorkspace gs_get_cytoframe
#' @importFrom flowWorkspace gs_get_pop_paths
#' @importFrom flowWorkspace gs_pop_add
#' @importFrom flowWorkspace gs_pop_get_count_fast
#' @importFrom flowWorkspace gs_pop_get_data
#' @importFrom flowWorkspace gs_pop_get_gate
#' @importFrom flowWorkspace gs_pop_get_parent
#' @importFrom flowWorkspace gs_pop_get_stats
#' @importFrom flowWorkspace gs_pop_remove
#' @importFrom flowWorkspace gs_pop_set_gate
#' @importFrom flowWorkspace gs_pop_set_name
#' @importFrom flowWorkspace load_cytoframe_from_fcs
#' @importFrom flowWorkspace load_cytoset_from_fcs
#' @importFrom flowWorkspace pop.MFI
#' @importFrom flowWorkspace recompute
#' @importFrom flowWorkspace transformerList
#' @importFrom fs path_home
#' @importFrom grDevices col2rgb
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices contourLines
#' @importFrom grDevices densCols
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom grDevices rgb
#' @importFrom graph graphNEL nodes
#' @importFrom htmlwidgets onRender
#' @importFrom magrittr "%>%"
#' @importFrom methods cbind2
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom methods rbind2
#' @importFrom rhandsontable hot_cols
#' @importFrom rhandsontable hot_context_menu
#' @importFrom rhandsontable hot_validate_numeric
#' @importFrom rhandsontable rHandsontableOutput
#' @importFrom rhandsontable renderRHandsontable
#' @importFrom rhandsontable rhandsontable
#' @importFrom scales alpha
#' @importFrom shinyFiles getVolumes
#' @importFrom shinyFiles parseDirPath
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles shinyDirButton
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyjs alert
#' @importFrom shinyjs delay
#' @importFrom shinyjs disable
#' @importFrom shinyjs enable
#' @importFrom shinyjs extendShinyjs
#' @importFrom shinyjs hide
#' @importFrom shinyjs inlineCSS
#' @importFrom shinyjs js
#' @importFrom shinyjs runjs
#' @importFrom shinyjs show
#' @importFrom shinyjs useShinyjs
#' @importFrom stats density
#' @importFrom stats rnorm
#' @importFrom stats setNames
#' @importFrom utils download.file
#' @importFrom utils tail
#' @importFrom utils write.csv
#' @importFrom utils write.table
#' @importMethodsFrom Rgraphviz AgNode
#' @importMethodsFrom flowWorkspace "["
#' @importMethodsFrom flowWorkspace "[["
#' @importMethodsFrom flowWorkspace "[[<-"
#' @importMethodsFrom flowWorkspace "colnames<-"
#' @importMethodsFrom flowWorkspace GatingSet
#' @importMethodsFrom flowWorkspace colnames
#' @importMethodsFrom flowWorkspace compensate
#' @importMethodsFrom flowWorkspace plot
#' @importMethodsFrom flowWorkspace transform

#Init Variables----
script <- "
$('#comp tr td').each(function() {
  var cellValue=$(this).text();
  if (cellValue == 100.0000000) {
    $(this).css('background-color', 'lightgrey');
    $(this).css('font-weight', 'normal')
  } else if (cellValue == 0) {
    $(this).css('font-weight', 'normal')
  } else if (cellValue > 0 && cellValue < 10) {
    $(this).css('background-color', '#FFFFB2');
    $(this).css('font-weight', 'normal')
  } else if (cellValue >= 10 && cellValue < 20) {
    $(this).css('background-color', '#FED976');
    $(this).css('font-weight', 'normal')
  } else if (cellValue >= 20 && cellValue < 40) {
    $(this).css('background-color', '#FEB24C');
    $(this).css('font-weight', 'normal')
  } else if (cellValue >= 40 && cellValue < 70) {
    $(this).css('background-color', '#FD8D3C');
    $(this).css('font-weight', 'normal')
  } else if (cellValue >= 70 && cellValue < 100) {
    $(this).css('background-color', '#F03B20');
    $(this).css('font-weight', 'normal')
  } else if (cellValue > 100) {
    $(this).css('background-color', '#BD0026');
    $(this).css('font-weight', 'normal')
  }
})"
changeHook <- "
function(el,x) {
  var hot=this.hot;
  var cellChanges=[];
  var changefn=function(changes,source) {
    if (source === 'edit'
      || source === 'undo'
      || source === 'autofill'
      || source === 'paste') {
      row=changes[0][0];
      col=changes[0][1];
      oldval=changes[0][2];
      newval=changes[0][3];
      if (oldval !== newval) {
        var cell=hot.getCell(row, col);
        cell.style.color='red';
        cellChanges.push({'rowid':row, 'colid':col});
      }
    }
  }
  var renderfn=function(isForced) {
    for(i=0; i < cellChanges.length; i++) {
      var rowIndex=cellChanges[i]['rowid'];
      var columnIndex=cellChanges[i]['colid'];
      var cell=hot.getCell(rowIndex, columnIndex);
      cell.style.color='red';
    }
  }
  var loadfn=function(initialLoad) {
    for(i=0; i < cellChanges.length; i++) {
      var rowIndex=cellChanges[i]['rowid'];
      var columnIndex=cellChanges[i]['colid'];
      var cell=hot.getCell(rowIndex, columnIndex);
      cell.style.color='black';
    }
    cellChanges=[] }
  hot.addHook('afterChange', changefn);
  hot.addHook('afterRender', renderfn);
  hot.addHook('afterLoadData', loadfn);
}"
jsCode <- "
shinyjs.disableTab=function(name) {
  var tab=$('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}
shinyjs.enableTab=function(name) {
  var tab=$('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}"
css <- "
.nav li a.disabled {
  color: #737373 !important;
    cursor: not-allowed !important;
}"
bc <- new.env()
tempEnv <- new.env()
bc$fluoCh <- ""
bc$Y <- ""
bc$prolifChannel <- NULL
bc$prolifParent <- NULL
bc$concat <- NULL
bc$gs <- NULL
bc$blackB <- "color:white; background-color:black; border-color:black"
bc$whiteB <- "color:black; background-color:white; border-color:black"
bc$sliderCol <- "+span>.irs-bar {background: black; border-top: black;
        border-bottom: black}"

#Internal funcs----
disableTabs <- function() {
    js$disableTab("ancestryTab")
    js$disableTab("overlayTab")
    js$disableTab("prolifTab")
    js$disableTab("tsneTab")
    js$disableTab("resultTab")
}

##File----
initializing <- function() {
    setwd(bc$filePath)
    cs <- load_cytoset_from_fcs(bc$fileList)
    bc$uncompGS <- GatingSet(cs)
    bc$gs <- gs_clone(bc$uncompGS)
    cf <- gs_get_cytoframe(bc$gs, bc$onlySampleIDs[1])
    bc$compDFs <- list(spillover(cf)[[1]])
    names(bc$compDFs) <- "Cytometer-defined"
    rownames(bc$compDFs[[1]]) <- colnames(bc$compDFs[[1]])
    compensate(bc$gs, bc$compDFs[[1]])
    colnames(bc$gs) <- colnames(cs)
    indices <- grepl("SC", colnames(bc$gs)) != TRUE
    preFluoChannels <- colnames(bc$gs)[indices]
    bc$fluoCh <- preFluoChannels[colnames(bc$gs)[indices] != "Time"]
    bc$customAxis <- rep(list(c(4.42, 0, -100)), length(bc$fluoCh))
    biexp <- flowjo_biexp_trans(channelRange=4096, maxValue=262144,
                                pos=bc$customAxis[[1]][1],
                                neg=bc$customAxis[[1]][2],
                                widthBasis=bc$customAxis[[1]][3])
    bc$gs <- transform(bc$gs, transformerList(bc$fluoCh, biexp))
    ff <- cytoframe_to_flowFrame(cs[[length(cs)]])
    bc$channelNames <- ff@parameters@data[,2]
    ch <- vapply(seq(colnames(bc$gs)), function(x)
        paste(colnames(bc$gs)[x], bc$channelNames[x], sep=" :: "),
        character(1))
    na <- vapply(bc$channelNames, function(x)
        is.na(x), TRUE)
    ch[na] <- colnames(bc$gs)[na]
    bc$ch <- as.list(colnames(bc$gs))
    names(bc$ch) <- ch
    completeFluoCh <- unlist(bc$ch)
    completeFluoCh <- completeFluoCh[completeFluoCh != "Time"]
    indices <- grepl("SC", completeFluoCh) != TRUE
    bc$completeFluoCh <- completeFluoCh[indices]
    bc$X <- colnames(bc$gs)[grep("FS", colnames(bc$gs))][1]
    bc$Y <- colnames(bc$gs)[grep("SS", colnames(bc$gs))][1]
    bc$sampleList <- bc$fileList[bc$onlySampleIDs]
    bc$results <- data.frame(matrix(ncol=1, nrow=length(bc$sampleList)))
    rownames(bc$results) <- substr(bc$sampleList, 1, nchar(bc$sampleList) - 4)
    colnames(bc$results) <- NA
    bc$hTable <- reactiveValues(d=as.data.frame(bc$compDFs[[1]]*100))
    bc$choices <- "Cytometer-defined"
    bc$appliedMatrix <- "Cytometer-defined"
}

initInpUpd <- function() {
    updateSelectInput(inputId="samp", choices=c("", bc$sampleList),
                      selected=bc$sampleList[1])
    updateSelectInput(inputId="previewSample", choices=c("", bc$fileList),
                      selected=bc$sampleList[1])
    updateSelectInput(inputId="bgPreviewSample", choices=c("", bc$sampleList))
    updateSelectInput(inputId="prolifSButt", choices=c("", bc$sampleList),
                      selected=bc$sampleList[1])
    updateSelectInput(inputId="Y", choices=c(bc$ch),
                      selected=bc$Y)
    updateSelectInput(inputId="X", choices=c(bc$ch),
                      selected=bc$X)
    updateSelectInput(inputId="ovY", choices=c("", bc$ch))
    updateSelectInput(inputId="ovX", choices=c("", bc$ch))
    updateSelectInput(inputId="rsCh", choices=c("", bc$completeFluoCh))
    js$enableTab("plotTab")
    js$enableTab("compTab")
}

##Axis----
addYAxis <- function(Y, axis, font) {
    index <- which(bc$fluoCh == Y)
    if(Y == "Time") {
        axisTicks(2, "time", bc$customAxis[[index]])
    } else if(grepl("SC", Y)) {
        axisTicks(2, "cont", bc$customAxis[[index]])
        if(axis == TRUE) {
            title(ylab=Y, cex.lab=(1+font/10), line=3.3, font.lab=2)
        }
    } else {
        axisTicks(2, "log", bc$customAxis[[index]])
        if(axis == TRUE) {
            id <- which(colnames(bc$gs) == Y)
            title(ylab=names(bc$ch)[id], cex.lab=(1+font/10), line=2.6,
                  font.lab=2)
        }
    }
}

axisTicks <- function(axisID, typ, customAxis) {
    if(typ == "time") {
        axisTTime(axisID)
    } else if(typ == "log") {
        axisTLog(axisID, customAxis)
    } else if(typ == "cont") {
        axisTCont(axisID)
    } else if(typ == "histo") {
        axisTHist(axisID)
    } else if(typ == "Overlaid histogram") {
        axisTOv(axisID)
    }
}

axisTTime <- function(axisID) {
    majorT <- axis(axisID, labels=FALSE, lwd=0)
    lastMajorInvT <- majorT[length(majorT)] + majorT[2] - majorT[1]
    minorT <- seq_len(lastMajorInvT/100)*100
    axis(axisID, at=majorT, lwd=2, labels=FALSE)
    axis(axisID, at=minorT, lwd=1, labels=FALSE, tck=-0.015)
    axis(axisID, at=majorT, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=majorT/100)
}

axisTLog <- function(axisID, customAxis) {
    trans <- flowjo_biexp(channelRange=4096, maxValue=262144,
                          pos=customAxis[1],
                          neg=customAxis[2],
                          widthBasis=customAxis[3])
    majorT <- trans(c(rev(10^seq_len(6)*-1), 0, 10^seq_len(6)))
    labels <- expression(10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1,
                         0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)
    neg <- majorT[seq_len(7)]
    pos <- majorT[7:length(majorT)]
    minorneg <- vapply(seq(pos), function(i)
        (pos[i+1]-pos[i])*log10(seq_len(9))+pos[i], numeric(9))
    minorpos <- vapply(seq(neg), function(i)
        (neg[i+1]-neg[i])*(1 - log10(rev(seq_len(9))))+neg[i], numeric(9))
    minorT <- list(as.list(minorneg[,-7]), as.list(minorpos[,-7]))
    minorT <- as.list(unlist(minorT))
    axis(axisID, at=minorT, lwd=2, labels=FALSE, tck=-0.015)
    adaptedLabels <- labels
    if(abs(majorT[6] - majorT[8]) < 600) {
        adaptedLabels[6] <- adaptedLabels[8] <- ""
    }
    if(abs(majorT[5] - majorT[9]) < 600) {
        adaptedLabels[5] <- adaptedLabels[9] <- ""
    }
    if(abs(majorT[4] - majorT[10]) < 600) {
        adaptedLabels[4] <- adaptedLabels[10] <- ""
    }
    if(abs(majorT[3] - majorT[11]) < 600) {
        adaptedLabels[3] <- adaptedLabels[11] <- ""
    }
    axis(axisID, at=majorT, lwd=2, labels=FALSE)
    axis(axisID, at=majorT, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=adaptedLabels)
}

axisTCont <- function(axisID) {
    majorT <- seq(0, 250000, by=50000)
    minorT <- seq(0, 260000, by=10000)
    axis(axisID, at=majorT, lwd=2, labels=FALSE)
    axis(axisID, at=minorT, lwd=1, labels=FALSE, tck=-0.015)
    axis(axisID, at=majorT, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=c("0", "50K", "100K", "150K", "200K", "250K"))
}

axisTHist <- function(axisID) {
    majorT <- axis(2, labels=FALSE, lwd=0)
    minorT <- seq(0, max(majorT)*1.2, by=majorT[1]+majorT[2]/5)
    axis(2, at=majorT, lwd=2, labels=FALSE)
    axis(2, at=minorT, lwd=1, labels=FALSE, tck=-0.015)
    shortL <- NULL
    for(i in seq(majorT)) {
        if(majorT[i] >= 1000) {
            n <- majorT[i]
            n2 <- n - as.numeric(substr(n, nchar(n)-2, nchar(n)))
            diminishedNumber <- n - n2
            if(diminishedNumber/1000 == 0) {
                shortL[[i]] <- paste0(n2/1000, "K")
            } else {
                shortL[[i]] <- paste0(n2/1000, ".", diminishedNumber/100, "K")
            }
        } else {
            shortL[[i]] <- majorT[i]
        }
    }
    axis(2, at=majorT, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=as.character(majorT))
}

axisTOv <- function(axisID) {
    majorT <- axis(2, labels=FALSE, lwd=0)
    axis(2, at=majorT, lwd=2, labels=FALSE)
    axis(2, at=seq(0, 1.2, by=majorT[1]+majorT[2]/4), lwd=1, tck=-0.015,
         labels=FALSE)
    axis(2, at=majorT, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=(majorT*100))
}

axisTtSNE <- function(xLim, yLim) {
    seqX1 <- seq(xLim[1]*0.9, xLim[2]*0.9, length.out=6)
    seqX2 <- seq(xLim[1]*0.9, xLim[2]*0.9, length.out=21)
    axis(1, at=seqX1, lwd=2, labels=FALSE)
    axis(1, at=seqX2, lwd=1, labels=FALSE, tck=-0.015)
    axis(1, at=seqX1, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=as.character(seq(0, 100, by=20)))
    seqY1 <- seq(yLim[1]*0.9, yLim[2]*0.9, length.out=6)
    seqY2 <- seq(yLim[1]*0.9, yLim[2]*0.9, length.out=21)
    axis(2, at=seqY1, lwd=2, labels=FALSE)
    axis(2, at=seqY2, lwd=1, labels=FALSE, tck=-0.015)
    axis(2, at=seqY1, lwd=0, cex.axis=1.1, line=-0.3, las=1,
         labels=as.character(seq(0, 100, by=20)))
    title(xlab="t-SNE 1", ylab="t-SNE 2", cex.lab=1.5, line=2.5, font.lab=2)
}

##Plot----
newPlot <- function(ID, X, Y, parent, typ, axis, font, bgDF=NULL, bg=FALSE) {
    ff <- cytoframe_to_flowFrame(gs_pop_get_data(bc$gs[ID], parent)[[1]])
    if(length(ff@exprs) > 0) {
        if(X == "Time") {
            xLim <- NULL
        } else if(grepl("SC", X)) {
            xLim <- c(0, 264000)
        } else {
            xLim <- c(0, 4100)
        }
        bc$dF <- as.data.frame.array(ff@exprs)
        dfX <- bc$dF[,X]
        dfX[dfX < 0] <- 0
        if(typ != "Histogram") {
            newDotPlot(dfX, X, Y, typ, axis, font, bgDF, xLim)
        } else {
            newHist(Y, dfX, axis, font, xLim)
        }
        index <- which(bc$fluoCh == X)
        if(X == "Time") {
            axisTicks(1, "time", bc$customAxis[[index]])
        } else if(grepl("SC", X)) {
            axisTicks(1, "cont", bc$customAxis[[index]])
        } else {
            axisTicks(1, "log", bc$customAxis[[index]])
        }
        if(axis == TRUE) {
            if(bg == TRUE) {
                lin <- 2
            } else {
                lin <- 3
            }
            id <- which(colnames(bc$gs) == X)
            title(xlab=names(bc$ch)[id], cex.lab=(1+font/10), line=lin,
                  font.lab=2)
        }
    } else {
        alert(paste0("The selected sample has 0 events under the selected
                     parent. Please choose another sample or change the
                     parent."))
    }
}

newDotPlot <- function(dfX, X, Y, typ, axis, font, bgDF, xLim) {
    dfY <- bc$dF[,Y]
    dfY[dfY < 0] <- 0
    if(Y == "Time") {
        yLim <- NULL
    } else if(grepl("SC", Y)) {
        yLim <- c(0, 264000)
    } else {
        yLim <- c(0, 4100)
    }
    if(typ == "Backgating") {
        plot(dfX, dfY, xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=0.5,
             lwd=0, col="grey", xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
        if(!is.null(bgDF)) {
            points(bgDF[,X], bgDF[,Y], col="red", pch=20, cex=0.5, lwd=0)
        }
    } else {
        if(typ != "Contour") {
            if(typ == "Density") {
                col <- c("grey0", "grey20", "grey40", "grey60", "grey80",
                         "grey100")
            } else {
                col <- c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
                         "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F",
                         "#9E0142")
            }
            colRamp <- colorRampPalette(col)
            density <- suppressWarnings(densCols(dfX, dfY, colramp=colRamp))
            plot(dfX, dfY, xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=0.5,
                 lwd=0, col=density, xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
        }
        if(typ == "Contour" || typ == "Pseudocolor + Contour") {

            addContour(dfX, dfY, X, Y, typ, xLim, yLim)
        }
    }
    addYAxis(Y, axis, font)
}

addContour <- function(dfX, dfY, X, Y, typ, xLim, yLim) {
    if(nrow(bc$dF) > 5000) {
        set.seed(1)
        sampledDF <- bc$dF[sample.int(nrow(bc$dF), 5000),]
    } else {
        sampledDF <- bc$dF
    }
    bc$z <- kde2d(sampledDF[,X], sampledDF[,Y], n=100)
    if(typ == "Contour") {
        plot(dfX, dfY, xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=0.5,
             lwd=0, col="black", xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
        for(i in contourLines(bc$z, nlevels=12)) {
            if(i$level != 5e-08) {
                polygon(i$x, i$y, col="white", lwd=1.2)
            }
        }
    } else {
        contour(bc$z, drawlabels=FALSE, xaxt="n", yaxt="n", add=TRUE,
                ann=FALSE, lwd=1.5, nlevels=12)
    }
}

newHist <- function(Y, dfX, axis, font, xLim) {
    if(length(dfX) > 1) {
        bc$referenceHist <- hist(dfX, breaks=20, plot=FALSE)
        refCounts <- bc$referenceHist$counts
        refDensity <- bc$referenceHist$density
        multiplier <- refCounts/refDensity
        bc$densLine <- stats::density(dfX)
        bc$densLine$y <- bc$densLine$y*multiplier[1]
        bc$maxDens <- max(bc$densLine[2]$y)
        if(Y == "Prolif") {
            multiplier <- 1.2
        } else {
            multiplier <- 1.05
        }
        plot(bc$densLine, xaxt="n", yaxt="n", ann=FALSE, xaxs="i", yaxs="i",
             xlim=xLim, ylim=c(0, bc$maxDens*multiplier))
        polygon(bc$densLine, col="grey", border="black", lwd=1)
        axisTicks(2, "histo", NULL)
    } else {
        plot(0, xaxt="n", yaxt="n", ann=FALSE, xaxs="i", yaxs="i", cex=0)
    }
    if(axis == TRUE) {
        splitRefCounts <- strsplit(as.character(max(refCounts)), "")[[1]]
        line <- length(splitRefCounts)*0.5 + 1.5
        title(ylab="Count", cex.lab=(1+font/10), line=line, font.lab=2)
    }
}

detectGate <- function(ID, X, Y, par, typ, obs, show, font, bg=FALSE) {
    popPaths <- gs_get_pop_paths(bc$gs)
    found <- NULL
    detectedTypes <- NULL
    if(length(popPaths) > 1) {
        for(i in seq(popPaths)[-1]) {
            gatesandPopsInfo <- gs_pop_get_gate(bc$gs, popPaths[i])[ID][[1]]
            detectedTypes[i] <- class(gatesandPopsInfo)[1]
            infoPar <- gatesandPopsInfo@parameters
            chX <- names(infoPar[1])
            if(typ != "Histogram") {
                if(length(infoPar) > 1) {
                    chY <- names(infoPar[2])
                } else {
                    chY <- "none"
                }
            } else {
                chY <- "Histogram" }
            if(X == chX && Y == chY
               || X == chY && Y == chX
               || typ == "Histogram" && length(infoPar) == 1 && X == chX) {
                detectedID <- which(popPaths == popPaths[i])
                preDetec <- gs_pop_get_count_fast(bc$gs[ID], "freq")
                detectedParent <- preDetec[detectedID - 1][,3][[1]]
                if(par == detectedParent) {
                    found[i-1] <- gatesandPopsInfo@filterId
                }
            }
        }
        if(obs == "plot" || obs == ">= 4") {
            bc$found <- found[!is.na(found)]
            detectedTypes <- detectedTypes[!is.na(detectedTypes)]
            if(obs == "plot" && !is.null(bc$found)) {
                plotGate(ID, X, Y, typ, detectedTypes, show,
                         font, bg)
            }
        }
    }
}

plotGate <- function(ID, X, Y, typ, types, show, font, bg=FALSE) {
    for(i in seq(bc$found)) {
        coords <- gs_pop_get_gate(bc$gs, bc$found[i])[[ID]]
        if(names(coords@parameters[1]) == X) {
            x <- 1
            y <- 2
        } else {
            x <- 2
            y <- 1
        }
        axisXLim <- par("usr")[2]
        axisYLim <- par("usr")[4]
        if(typ != "Histogram") {
            gateHist(coords, x, y, axisXLim, show)
            if(is(coords,"rectangleGate")) {
                rect(xleft=bc$xLeft, xright=bc$xRight,
                     ybottom=bc$yBottom, ytop=bc$yTop, lwd=3)
            }
        } else {
            segments(coords@min[[x]], axisYLim*0.75, coords@max[[x]],
                     axisYLim*0.75, lwd=3)
            segments(coords@min[[x]], axisYLim*0.72, coords@min[[x]],
                     axisYLim*0.78, lwd=3)
            segments(coords@max[[x]], axisYLim*0.72, coords@max[[x]],
                     axisYLim*0.78, lwd=3)
        }
        popPaths <- gs_get_pop_paths(bc$gs)
        popShortPaths <- gs_get_pop_paths(bc$gs, path=1)
        subpop <- popPaths[which(popShortPaths == bc$found[i])]
        gatedPopInfo <- gs_pop_get_count_fast(bc$gs[ID],
                                              subpopulations=subpop)
        preFreq <- gatedPopInfo[,4][[1]]*100/gatedPopInfo[,5][[1]]
        freq <- paste0(sprintf("%.1f", preFreq), "%")
        gateLeg(ID, X, Y, typ, show, font, bg, freq, coords, x, i)
    }
}

gateHist <- function(coords, x, y, axisXLim, show) {
    if(is(coords, "polygonGate")) {
        bound <- coords@boundaries
        bc$xLeft <- min(bound[,1])
        bc$xRight <- max(bound[,1])
        bc$yBottom <- min(bound[,2])
        bc$yTop <- max(bound[,2])
        for(j in seq(nrow(bound))) {
            if(j > 1) {
                segments(bound[,1][j-1], bound[,2][j-1], bound[,1][j],
                         bound[,2][j], lwd=3)
            }
        }
    } else {
        if(coords@min[[x]] == "-Inf") {
            bc$xLeft <- -1000
        } else {
            bc$xLeft <- coords@min[[x]]
        }
        if(coords@max[[x]] == "Inf") {
            bc$xRight <- par("usr")[2]*1.1
        } else {
            bc$xRight <- coords@max[[x]]
        }
        if(coords@min[[y]] == "-Inf") {
            bc$yBottom <- -1000
        } else {
            bc$yBottom <- coords@min[[y]]
        }
        if(coords@max[[y]] == "Inf") {
            axisYLim <- par("usr")[4]
            bc$yTop <- axisYLim*1.1
        } else {
            bc$yTop <- coords@max[[y]]
        }
    }
    if((bc$yBottom+bc$yTop)/2 < max(axis(y, labels=FALSE, lwd=0))/2) {
        if(show == TRUE) {
            bc$yJust <- -0.5
        } else {
            bc$yJust <- -1.2
        }
    } else {
        if(show == TRUE) {
            bc$yJust <- 1.5
        } else {
            bc$yJust <- 2.2
        }
    }
}

gateLeg <- function(ID, X, Y, typ, show, font, bg, freq, coords, x, i) {
    if(show == TRUE) {
        title <- bc$found[i]
    } else {
        title <- NULL
    }
    col <- rgb(1,1,1, 0.6)
    axisYLim <- par("usr")[4]
    if(length(bc$found) == 4) {
        if(bg == TRUE) {
            qXCoords <- c(0.2, 0.8, 0.8, 0.2)
            qYCoords <- c(1.1, 1.1, 0.5, 0.5)
        } else {
            if(names(coords@min)[1] == X) {
                qXCoords <- c(0.1, 0.9, 0.9, 0.1)
                qYCoords <- c(1.05, 1.05, 0.25, 0.25)
            } else {
                qXCoords <- c(0.9, 0.9, 0.1, 0.1)
                qYCoords <- c(0.25, 1.05, 1.05, 0.25)
            }
        }
        if(show == TRUE) {
            if(bc$found[i] == paste0("Q1: ", X, "- ", Y, "+")
               || bc$found[i] == paste0("Q2: ", X, "+ ", Y, "+")
               || bc$found[i] == paste0("Q3: ", X, "+ ", Y, "-")
               || bc$found[i] == paste0("Q4: ", X, "- ", Y, "-")) {
                title <- substring(bc$found[i], 1, 2)
            }
        }
        legX <- par("usr")[2]*qXCoords[i]
        legY <- axisYLim*qYCoords[i]
        legend(legX, legY, title=title, legend=freq, cex=1+font/10, bg=col,
               box.lwd=0, x.intersp=-0.5, y.intersp=0.8, text.font=2,
               xjust=0.5, yjust=1.5)
    } else {
        if(typ != "Histogram") {
            legX <- (bc$xLeft+bc$xRight)/2
            legY <- (bc$yBottom+bc$yTop)/2
            legend(legX, legY, title=title, legend=freq, cex=1+font/10, bg=col,
                   box.lwd=0, x.intersp=-0.5, y.intersp=0.8, text.font=2,
                   xjust=0.5, yjust=bc$yJust)
        } else {
            legX <- (coords@min[[x]]+coords@max[[x]])/2
            legY <- axisYLim*0.95
            legend(legX, legY, title=title, legend=freq, cex=1+font/10, bg=col,
                   box.lwd=0, x.intersp=-0.5, y.intersp=0.8, text.font=2,
                   xjust=0.5)
        }
    }
}

hideTools <- function() {
    hide("rectang", TRUE, "fade")
    hide("polygon", TRUE, "fade")
    hide("quad", TRUE, "fade")
    hide("interval", TRUE, "fade")
}

secStep <- function() {
    delay(500, shinyjs::show("gateOk", TRUE, "fade"))
    delay(500, shinyjs::show("gateCancel", TRUE, "fade"))
    delay(500, shinyjs::show("drawGateName", TRUE, "fade"))
    disable("gateOk")
}

reAddPops <- function() {
    setwd(bc$filePath)
    popPaths <- gs_get_pop_paths(bc$gs)
    tempEnv$gs <- gs_clone(bc$uncompGS)
    compensate(tempEnv$gs, bc$compDFs[bc$appliedMatrix][[1]])
    popParents <- gs_pop_get_count_fast(bc$gs[1], "freq")[,3][[1]]
    popGates <- vapply(popPaths[-1], function(x)
        list(gs_pop_get_gate(bc$gs, x)[1][[1]]),
        list(length(popPaths[-1])))
    compensate(tempEnv$gs, bc$compDFs[bc$appliedMatrix][[1]])
    for(i in seq(popGates)) {
        gs_pop_add(tempEnv$gs, popGates[[i]], parent=popParents[i])
    }
    recompute(tempEnv$gs)
}

##Compensation----
compGen <- function(ID, newMatrix, ov, parent="root") {
    setwd(bc$filePath)
    plotN <- length(bc$fluoCh) - 1
    mfRow <- c(plotN, plotN)
    par(oma=c(1, 2, 2, 1), mar=c(0, 0, 0, 0), mfrow=mfRow, lwd=2)
    cf <- load_cytoframe_from_fcs(bc$fileList[ID])
    originalMatrix <- bc$compDFs[[1]]
    rownames(originalMatrix) <- colnames(originalMatrix)
    if(is.null(newMatrix)) {
        compensate(cf, originalMatrix)
    } else {
        compensate(cf, newMatrix)
    }
    for(i in seq(bc$fluoCh)) {
        trans <- flowjo_biexp_trans(channelRange=4096, maxValue=262144,
                                    pos=bc$customAxis[[i]][1],
                                    neg=bc$customAxis[[i]][2],
                                    widthBasis=bc$customAxis[[i]][3])
        cf <- transform(cf, transformerList(bc$fluoCh[i], trans))
    }
    df <- as.data.frame(cytoframe_to_flowFrame(cf)@exprs)
    uncDF <- as.data.frame.array(cytoframe_to_flowFrame(
        gs_pop_get_data(bc$uncompGSfortable[ID], parent)[[1]])@exprs)
    col <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
             "#FC4E2A", "#E31A1C", "#BD0026", "#800026")
    col <- alpha(col, alpha=0.4)
    byRows <- vapply(seq(plotN), function(i)
        list(as.numeric(newMatrix[i,])[seq_len(i)*-1]*100), list(plotN))
    byCols <- vapply(seq(plotN), function(i)
        list(as.numeric(newMatrix[,i])[seq_len(i)*-1]*100), list(plotN))
    same <- vapply(seq(unlist(byRows)), function(i)
        list(c(unlist(byRows)[i], unlist(byCols)[i])), list(2))
    range <- vapply(seq_len(11), function(i)
        list(c(9*i-9, 9*i)), list(2))
    lim <- c(0, 4100)
    id <- 0
    for(i in seq(plotN)) {
        for(j in seq(plotN + 1 - i)) {
            id <- id + 1
            compPlot(i, j, same, range, lim, id, plotN,
                     df, uncDF, col, ov)
        }
        if(i != plotN) {
            for(k in seq(i)) {
                plot.new()
            }
        }
    }
}

compPlot <- function(i, j, same, range, lim, id, plotN, df, uncDF, col, ov) {
    X <- bc$fluoCh[j+i]
    Y <- bc$fluoCh[i]
    set.seed(3)
    sampledDF <- df[sample.int(nrow(df), 1000),]
    uncsampledDF <- uncDF[sample.int(nrow(df), 1000),]
    dfX <- sampledDF[,X]
    dfY <- sampledDF[,Y]
    plot(dfX, dfY, xaxt="n", yaxt="n", ann=FALSE, pch=20,
         cex=0.5, lwd=0, xlim=lim, ylim=lim, xaxs="i", yaxs="i")
    prod <- same[[id]]
    for(l in seq_len(11)) {
        if(max(prod) == 0) {
            bgCol <- "white"
        } else {
            if(max(abs(prod)) > range[[l]][1]
               && max(abs(prod)) <= range[[l]][2]
               && max(abs(prod)) <= 81) {
                bgCol <- col[l]
            } else if((max(abs(prod)) > 81)) {
                bgCol <- col[9]
            }
        }
    }
    rLim <- par("usr")
    rect(rLim[1], rLim[3], rLim[2], rLim[4], col=bgCol)
    points(dfX, dfY, xaxt="n", yaxt="n", ann=FALSE, pch=20,
           cex=0.5, lwd=0, xlim=lim, ylim=lim, xaxs="i", yaxs="i")
    if(ov == TRUE) {
        points(uncsampledDF[,X], uncsampledDF[,Y], col="grey",ann=FALSE,
               xaxt="n", yaxt="n", pch=20, cex=0.5, lwd=0,
               xlim=lim, ylim=lim, xaxs="i", yaxs="i")
    }
    if(length(bc$fluoCh) >= 10) {
        textSize <- 0.6
    } else if(length(bc$fluoCh) < 10 && length(bc$fluoCh) > 4) {
        textSize <- rev(5:9*0.12)[length(bc$fluoCh) - 4]
    } else if(length(bc$fluoCh) <= 4) {
        textSize <- 1.2
    }
    if(i == 1) {
        mtext(X, side=3, padj=-0.5, cex=textSize)
    }
    if(j == 1 || i > j && i == j - 1) {
        mtext(Y, side=2, padj=-0.5, cex=textSize)
    }
}

matChoices <- function() {
    bc$choices <- names(bc$compDFs)
    names(bc$choices) <- bc$choices
    allNames <- c()
    for(i in seq(bc$choices)) {
        if(bc$choices[i] == bc$appliedMatrix) {
            allNames[i] <- paste(bc$appliedMatrix, "(applied)")
        } else {
            allNames[i] <- bc$choices[i]
        }
    }
    names(bc$choices) <- allNames
}

##Ancestry----
ancestryGen <- function(typ, previewSamp, bgPop) {
    par(mfrow=c(3, 5), mar=c(4,5,2,1) + 0.1, lwd=2)
    if(typ != "" && previewSamp != "") {
        ID <- match(previewSamp, bc$fileList)
        pop <- bgPop
        if(typ == "Backgating" && pop != "") {
            cs <- gs_pop_get_data(bc$gs[ID], pop)
            ff<- cytoframe_to_flowFrame(cs[[1]])
            bgDF <- as.data.frame.array(ff@exprs)
        } else {
            bgDF <- NULL
        }
        popPaths <- gs_get_pop_paths(bc$gs)
        gates <- vapply(popPaths[-1], function(x)
            list(gs_pop_get_gate(bc$gs, x)[[1]]), list(1))
        gateChs <- vapply(gates, function(x)
            list(names(x@parameters)), list(1))
        seq1 <- seq(gateChs)[duplicated(gateChs)]
        seq2 <- seq(gateChs)[duplicated(gateChs, fromLast=TRUE)]
        duplicateIDs <- sort(unique(c(seq1, seq2)))
        data <- gs_pop_get_count_fast(bc$gs[[1]])
        duplicateParents <- data[,"Parent"][[1]][duplicateIDs]
        if(length(duplicateIDs) > 0) {
            index <- duplicateIDs[duplicated(duplicateParents)]
            plotSeq <- seq(gateChs)[-index]
        } else {
            plotSeq <- seq(gateChs)
        }
        dat <- data.frame(x=numeric(0), y=numeric(0))
        acPlotGen(dat, plotSeq, gates, popPaths, typ, ID, bgDF, pop)
        enable("exportImageAncestry")
    } else {
        disable("bgPop")
    }
}

acPlotGen <- function(dat, plotSeq, gates, popPaths, typ, ID, bgDF, pop) {
    max <- length(plotSeq)/10
    withProgress(message="Generating", detail="plot 0", value=0, max=max, {
        for(i in (plotSeq+1)) {
            dat <- rbind(dat, data.frame(x=rnorm(1), y=rnorm(1)))
            detail <- paste("plot", which((plotSeq+1) == i))
            incProgress(0.1, detail=detail)
            Sys.sleep(0.1)
            gateInfo <- gates[popPaths[i]][[1]]
            xCh <- names(gateInfo@parameters[1])
            if(length(gateInfo@parameters) > 1) {
                yCh <- names(gateInfo@parameters[2])
                type <- typ
            } else {
                yCh <- "Histogram"
                type <- "Histogram"
            }
            shortPath <- gs_get_pop_paths(bc$gs, path=1)
            index <- which(shortPath == gateInfo@filterId)
            gateFullName <- popPaths[index]
            if(i > 2) {
                stSp <- strsplit(gateFullName, split="/")
                end <- tail(stSp[[1]], 2)[1]
                i2 <- which(shortPath == end)
                parent <- popPaths[i2]
            } else {
                parent <- "root"
            }
            newPlot(ID, xCh, yCh, parent, type, TRUE, 2, bgDF, TRUE)
            detectGate(ID, xCh, yCh, parent, type, "plot", TRUE, 2, TRUE)
            if(parent == "root") {
                parent <- "ungated"
            } else {
                parent <- shortPath[which(popPaths == parent)]
            }
            if(typ == "Backgating" && pop != "") {
                if(!is.null(bgDF)) {
                    if(parent == pop) {
                        title(main=list(parent, cex=1.5, col="red"))
                    } else {
                        title(main=list(parent, cex=1.5, col="black"))
                    }
                }
            } else {
                title(main=list(parent, cex=1.5, col="black"))
            }
        }
    })
}

##Overlays----
overlay <- function(IDs, X, Y, parent, typ, tone, axis, font) {
    cfs <- vapply(IDs, function(i)
        list(gs_pop_get_data(bc$gs[i], parent)[[1]]), list(length(IDs)))
    shortName <- substr(bc$fileList, 1, nchar(bc$fileList) - 4)
    bc$sampleNam <- shortName[IDs]
    if(all(sapply(cfs, function(x) nrow(x) > 0))) {
        par(mar=c(4, 6, 1, 25), lwd=2)
        if(X == "Time") {
            bc$xLim <- NULL
        } else if(grepl("SC", X)) {
            bc$xLim <- c(0, 264000)
        } else {
            bc$xLim <- c(0, 4100)
        }
        setColor(IDs, typ, tone)
        dfs <- vapply(cfs, function(x)
            list(as.data.frame(cytoframe_to_flowFrame(x)@exprs)),
            list(length(cfs)))
        dfsX <- vapply(dfs, function(x)
            list(x[,X]), list(length(dfs)))
        dfsX <- lapply(dfsX, function(x) {
            x[x < 0] <- 0
            x
        })
        if(typ == "Dot plot") {
            ovDotPlot(IDs, Y, dfs, dfsX)
        } else {
            ovHist(IDs, dfsX, typ)
        }
        multip <- 1.2
        ovAxis(IDs, X, Y, typ, axis, font, multip)
    } else {
        alert(paste0("The selected sample has 0 events under the selected
                     parent. Please choose another sample or change the
                     parent."))
    }
}

setColor <- function(IDs, typ, tone) {
    if(typ == "Overlaid histogram") {
        transparency <- 0.6
        bc$colorList <- c("grey"="transparent", "transparent"="black",
                          "transparent"="black")[seq_len(length(IDs))]
        bc$lty <- c(0, 1, 2)
        bc$lwd <- c(1, 3, 2)
    } else if(typ == "Offset histogram") {
        transparency <- 0.5
        bc$yLim <- c(0, 0)
        bc$offsets <- c(0)
        if(tone == "Colorful") {
            preColors <- c("grey", "red2", "green4", "blue3", "orange",
                           "black", "wheat3", "brown", "orchid3", "salmon")
        } else {
            prepreColors <- c("grey0", "grey10", "grey20", "grey30", "grey40",
                              "grey50", "grey60", "grey70", "grey80",
                              "grey90")
            multip <- length(prepreColors)/length(IDs)
            preColors <- vapply(seq(IDs), function(i)
                prepreColors[multip*i], character(1))
            preColors[1] <- prepreColors[1]
            preColors[length(IDs)] <- prepreColors[length(prepreColors)]
        }
        alphaColors <- vapply(seq(IDs), function(i)
            rgb(col2rgb(preColors[i])[1],
                col2rgb(preColors[i])[2],
                col2rgb(preColors[i])[3],
                255*transparency, maxColorValue=255),
            character(1))
        bc$colorList <- rep("black", length(IDs))
        names(bc$colorList) <- alphaColors
        bc$lty <- bc$lwd <- rep(1, length(IDs))
    } else if (typ == "Dot plot"){
        if(tone == "Colorful") {
            bc$colorList <- c("black", "red2")
        } else {
            bc$colorList <- c("black", "grey")
        }
    }
}

ovDotPlot <- function(IDs, Y, dfs, dfsX) {
    dfsY <- vapply(dfs, function(y)
        list(y[,Y]), list(length(dfs)))
    dfsY <- lapply(dfsY, function(y){
        y[y < 0] <- 0
        y
    })
    if(Y == "Time") {
        bc$yLim <- NULL
    } else if(grepl("SC", Y)) {
        bc$yLim <- c(0, 264000)
    } else {
        bc$yLim <- c(0, 4100)
    }
    dfFinal <- data.frame("X"=dfsX[[1]], "Y"=dfsY[[1]],
                          "color"=bc$colorList[1])
    dfFinal <- rbind(dfFinal, data.frame("X"=dfsX[[2]], "Y"=dfsY[[2]],
                                         "color"=bc$colorList[2]))
    set.seed(2)
    shuffledDF <- dfFinal[sample(nrow(dfFinal)), ]
    plot(shuffledDF[[1]], shuffledDF[[2]], xaxt="n", yaxt="n",
         ann=FALSE, pch=20, cex=0.5, lwd=0, col=shuffledDF[[3]],
         xlim=bc$xLim, ylim=bc$yLim, xaxs="i", yaxs="i")
}

ovHist <- function(IDs, dfsX, typ) {
    for(i in seq(IDs)) {
        densLine <- density(dfsX[[i]])
        densLine$y <- densLine$y/max(densLine$y)
        maxDens <- max(densLine$y)
        if(typ == "Overlaid histogram") {
            bc$yLim <- c(0, maxDens*1.05)
            if(i == 1) {
                plot(densLine, col="transparent", xaxt="n", yaxt="n",
                     ann=FALSE, xlim=bc$xLim, ylim=c(0, maxDens*1.05),
                     xaxs="i", yaxs="i")
            }
            polygon(densLine, col=names(bc$colorList[i]),
                    border=bc$colorList[i], lwd=bc$lwd[i], lty=bc$lty[i])
        } else {
            if(i != 1) {
                bc$offsets[i] <- bc$yLim[2]
            }
            bc$yLim[2] <- bc$yLim[2] + maxDens*0.7
        }
    }
    if(typ == "Offset histogram") {
        bc$yLim[2] <- bc$yLim[2] + maxDens*0.7/2*1.1
        plot(0, xaxt="n", yaxt="n", ann=FALSE, xlim=bc$xLim, ylim=bc$yLim,
             xaxs="i", yaxs="i", cex=0)
        for(i in rev(seq(IDs))) {
            densLine <- stats::density(dfsX[[i]])
            densLine$y <- densLine$y/max(densLine$y) + bc$offsets[i]
            polygon(densLine, col=names(bc$colorList[i]), border=NA, lwd=1.5,
                    lty=1)
            lines(densLine, lwd=1.5, lty=1)
        }
    }
}

ovAxis <- function(IDs, X, Y, typ, axis, font, multip) {
    if(X == "Time") {
        axisTicks(1, "time", NULL)
    } else if(grepl("SC", X)) {
        axisTicks(1, "cont", NULL)
    } else {
        axisTicks(1, "log", bc$customAxis[[which(bc$fluoCh == X)]])
    }
    if(typ == "Dot plot") {
        if(Y == "Time") {
            distance <- 2.8
            axisTicks(2, "time", NULL)
        } else if(grepl("SC", Y)) {
            distance <- 3.3
            axisTicks(2, "cont", NULL)
        } else {
            distance <- 2.8
            axisTicks(2, "log", bc$customAxis[[which(bc$fluoCh == Y)]])
        }
    } else if(typ == "Overlaid histogram") {
        axisTicks(2, typ, NULL)
        pch <- c(22, NA, NA)[seq_len(length(IDs))]
        lty <- c(0, 1, 2)[seq_len(length(IDs))]
    }
    if(axis == TRUE) {
        xLab <- names(bc$ch)[which(colnames(bc$gs) == X)]
        title(xlab=xLab, cex.lab=(1+font/10), line=3, font.lab=2)
        if(typ == "Overlaid histogram") {
            title(ylab="% of max", cex.lab=(1+font/10), line=2.5, font.lab=2)
        } else if(typ != "Offset histogram") {
            yLab <- names(bc$ch)[which(colnames(bc$gs) == Y)]
            title(ylab=yLab, cex.lab=(1+font/10), line=distance, font.lab=2)
        }
    }
    if(typ != "Overlaid histogram") {
        pch <- 22
        lty <- 0
    }
    if(typ == "Overlaid histogram" || typ == "Offset histogram") {
        ptBg <- names(rev(bc$colorList))
    } else {
        ptBg <- rev(bc$colorList)
    }
    leg <- legend("topright", bty="n", inset=c(-0.25, -0.025), pch=rev(pch),
                  text.width=strwidth("100"), lty=rev(lty), lwd=2, cex=1.3,
                  pt.cex=4, y.intersp=1.6, pt.bg=ptBg, col=rev(bc$colorList),
                  legend=rep(NA, length(bc$sampleNam)), xpd=TRUE)
    text(strwidth(rev(bc$sampleNam))*multip+leg$text$x*1.01, leg$text$y,
         rev(bc$sampleNam), pos=2, xpd=NA, font=2, cex=multip)
}

##Proliferation----
prolifGen <- function(prolifSB, show, font, ready, lab, grid) {
    undividedGate <- gs_pop_get_gate(bc$gs, "undivided")[[1]]
    bc$prolifChannel <- names(undividedGate@parameters)
    popPaths <- gs_get_pop_paths(bc$gs)
    popShortPaths <- gs_get_pop_paths(bc$gs, path=1)
    fullUndivided <- popPaths[which(popShortPaths == "undivided")]
    bc$prolifParent <- gs_pop_get_count_fast(bc$gs[bc$sampleList[1]])
    ID <- which(bc$prolifParent[,2] == fullUndivided)
    bc$prolifParent <- bc$prolifParent[,3][ID][[1]]
    par(mar=c(4,6,1,1) + 0.1, lwd=2)
    ID <- match(prolifSB, bc$fileList)
    newPlot(ID, bc$prolifChannel, "Prolif", bc$prolifParent, "Histogram",
            show, font)
    shinyjs::show("prolifSButt")
    shinyjs::show("prevProlifS")
    shinyjs::show("nextProlifS")
    if(ready == FALSE) {
        prolifNotReady(undividedGate)
    } else {
        hide("applyProlif")
        disable("step1")
        disable("step2")
        shinyjs::show("exportTableProlif")
        shinyjs::show("exportImageProlif")
        shinyjs::show("prolifLabel")
        shinyjs::show("prolifGrid")
        shinyjs::show("prolifTable")
        enable("step3")
        dfX <- bc$dF[,bc$prolifChannel]
        ref <- hist(dfX, breaks=150, plot=FALSE)
        prolifReady(dfX, ref, lab, grid)
    }
}

prolifNotReady <- function(undividedGate) {
    prolifOff()
    disable("step1")
    shinyjs::show("prolifSButt")
    shinyjs::show("prevProlifS")
    shinyjs::show("nextProlifS")
    shinyjs::show("applyProlif")
    enable("step2")
    undividedLimits <- c(undividedGate@min[[1]], undividedGate@max[[1]])
    ID <- bc$referenceHist$breaks > undividedLimits[1]
    approxD0X <- bc$referenceHist$breaks[ID]
    approxD0X <- approxD0X[approxD0X < undividedLimits[2]]
    peakID <- which(diff(sign(diff(bc$densLine$y))) == -2)
    xPeaks <- bc$densLine$x[peakID]
    yPeaks <- bc$densLine$y[peakID]
    ID <- which(diff(sign(diff(bc$densLine$y))) == 2)
    bc$betweenPeaks <- bc$densLine$x[ID]
    ID <- bc$betweenPeaks < mean(approxD0X)*1.05
    bc$betweenPeaks <- bc$betweenPeaks[ID]
    bc$refPeaks <- data.frame(x=NA, y=NA)
    for(i in seq(xPeaks)) {
        bc$refPeaks[i,] <- c(xPeaks[[i]], yPeaks[[i]])
    }
    IDs <- which(bc$refPeaks[[1]] > mean(approxD0X)*1.05)
    if(length(IDs) > 0) {
        bc$refPeaks <- bc$refPeaks[-IDs,]
    }
    bc$adjustedL <- c()
    for(i in seq(nrow(bc$refPeaks))) {
        if(i == 1) {
            mean <- mean(bc$refPeaks[,1]-bc$refPeaks[,1][1])/3
            bc$adjustedL[i] <- bc$betweenPeaks[1] - mean
        } else if(i == nrow(bc$refPeaks)) {
            bc$adjustedL[i] <- bc$refPeaks[,1][i]
        } else {
            mean <- mean(c(bc$betweenPeaks[i], bc$betweenPeaks[i-1]))
            bc$adjustedL[i] <- mean
        }
    }
    for(i in seq(bc$refPeaks[,1])) {
        text(rev(bc$adjustedL)[i], bc$maxDens*1.1, labels=i-1, cex=1.5)
        abline(v=bc$betweenPeaks[i], lty=2)
    }
}

prolifReady <- function(dfX, ref, lab, grid) {
    gap <- mean(bc$refPeaks[,1]-bc$refPeaks[,1][1])/nrow(bc$refPeaks)
    bc$gatecoords <- list()
    for(i in seq(bc$refPeaks[,1])) {
        if(lab == FALSE) {
            dX <- c(bc$refPeaks[,1][i]-gap, bc$refPeaks[,1][i]+gap)
            dBreaks <- ref$breaks[ref$breaks > dX[1]]
            dBreaks <- dBreaks[dBreaks < dX[2]]
            dCounts <- c()
            for(j in seq(dBreaks)) {
                ID <- which(ref$breaks == dBreaks[j])
                dCounts[j] <- ref$counts[ID]
            }
            y <- mean(dCounts)*13
            if(!is.na(y)) {
                if(y > bc$maxDens*1.1) {
                    y <- bc$maxDens*1.1
                }
            }
            labels <- abs(i-nrow(bc$refPeaks))
            text(bc$adjustedL[i], y, labels=labels, cex=1.5)
        } else {
            y <- bc$maxDens*1.1
            text(rev(bc$adjustedL)[i], y, labels=i-1, cex=1.5)
        }
        if(grid == TRUE) {
            abline(v=bc$betweenPeaks[i], lty=2)
        }
        if(i == 1) {
            bc$gatecoords[[i]] <- c(0, bc$betweenPeaks[1])
        } else if(i == nrow(bc$refPeaks)) {
            bc$gatecoords[[i]] <- c(bc$betweenPeaks[i-1], 4100)
        } else {
            bc$gatecoords[[i]] <- c(bc$betweenPeaks[i-1],
                                    bc$betweenPeaks[i])
        }
    }
}

prolifTableGen <- function(prolifSB) {
    ID <- match(prolifSB, bc$fileList)
    ff <- cytoframe_to_flowFrame(gs_pop_get_data(bc$gs[ID],
                                                 bc$prolifParent)[[1]])
    dF <- as.data.frame.array(ff@exprs)
    dfX <- dF[,bc$prolifChannel]
    bc$divPercents <- c()
    bc$divCounts <- c()
    bc$total <- length(dfX)
    for(i in seq(bc$gatecoords)) {
        filter <- dfX[dfX > bc$gatecoords[[i]][1]]
        filter <- filter[filter < bc$gatecoords[[i]][2]]
        bc$divCounts[i] <- length(filter)
        bc$divPercents[i] <- length(filter)*100/bc$total
    }
    bc$divCounts <- rev(bc$divCounts)
    round <- round(bc$divPercents, 2)
    bc$divPercents <- rev(as.numeric(format(round, nsmall=2)))
    bc$totalDivCount <- bc$total-bc$divCounts[1]
    bc$totalDivPerc <- 100-bc$divPercents[1]
    l0 <- "Proliferation assay data:"
    l1 <- paste("Number of divisions:", (nrow(bc$refPeaks) - 1))
    l2 <- paste("Divided cells:", bc$totalDivCount)
    l3 <- paste("Undivided cells:", bc$divCounts[1])
    lines <- c(l1, l2, l3)
    for(i in 2:nrow(bc$refPeaks)) {
        lines[i+2] <- paste0("Division ", (i - 1), ": ", bc$divCounts[i])
    }
    number <- as.numeric(format(round(bc$totalDivPerc, 2), nsmall=2))
    lines[length(lines) + 1] <- paste("Percent divided:", number)
    number <- bc$divPercents[1]
    lines[length(lines) + 1] <- paste("Percent undivided:", number)
    constantlinesln <- length(lines)
    for(i in 2:nrow(bc$refPeaks)) {
        lines[i+(constantlinesln-1)] <- paste0("Percent division ", (i - 1),
                                               ": ", bc$divPercents[i])
    }
    part2 <- paste(lines, collapse="<br/>", sep="<br/>")
    HTML("<br/>", paste(strong(l0)), "<br/>", part2)
}

prolifOff <- function() {
    hide("exportTableProlif")
    hide("exportImageProlif")
    hide("prolifLabel")
    hide("prolifGrid")
    hide("prolifTable")
    disable("step3")
}

##t-SNE----
tSNEGen <- function(dots, mode, pops, gOrS, ID, tSHighl, react, tIDs, tPops) {
    par(mar=c(4, 6, 3, 26), lwd=2)
    currentEnttS <- bc$entiretSNE
    x <- currentEnttS[,1]
    y <- currentEnttS[,2]
    xLim <- c(min(x)*1.1, max(x)*1.1)
    yLim <- c(min(y)*1.1, max(y)*1.1)
    nrows <- nrow(currentEnttS)
    cex <- dots/10
    if(mode == "Overlay Groups or Samples") {
        if(gOrS == "All") {
            enable("savetSNEPlot")
            plot(x, y, col="black", xaxt="n", yaxt="n", ann=FALSE, pch=20,
                 cex=cex, lwd=0, xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
            axisTtSNE(xLim, yLim)
        } else {
            if(!is.null(react)) {
                bc$tIDs <- tIDs
                tSGenGroupORS(currentEnttS, nrows, cex, xLim, yLim,
                              as.integer(react), gOrS)
            }
        }
        hide("showingtSNEEvents")
    } else {
        if(gOrS == "All") {
            dfRows <- c(1, nrows)
            currenttS <- currentEnttS
        } else {
            if(gOrS == "Group") {
                interval <- nrows/length(bc$concatSamples)
            } else if(gOrS == "Sample") {
                interval <- nrows/length(unlist(bc$concatSamples))
            }
            num <- interval*as.numeric(ID)
            dfRows <- c(num - interval + 1, num)
            currenttS <- currentEnttS[dfRows[1]:dfRows[2],]
        }
        if(mode == "Heatmap") {
            tSGenHeat(cex, xLim, yLim, tSHighl, dfRows, currenttS)
            shinyjs::show("showingtSNEEvents")
        } else if(mode == "Overlay Populations"
                  && !is.null(pops)) {
            tSGenOvPo(dfRows, currentEnttS, cex, xLim, yLim, tPops)
        }
        bc$toFormat <- as.integer((dfRows[2] - dfRows[1]) + 1)
    }
}

tSGenGroupORS <- function(currentEnttS, nrows, cex, xLim, yLim, ids, gOrS) {
    if(gOrS == "Group") {
        interval <- nrows/length(bc$concatSamples)
        listGroupORSamp <- names(bc$tSNEListofGroups[ids])
    } else if(gOrS == "Sample") {
        interval <- nrows/length(unlist(bc$concatSamples))
        listGroupORSamp <- names(bc$tSNEListofSamples[ids])
    }
    firs <- interval*as.numeric(bc$tIDs[1])
    sec <- interval*as.numeric(bc$tIDs[2])
    dfRows <- list(c(firs - interval + 1, firs),
                   c(sec - interval + 1, sec)
    )
    currenttS <- list(
        currentEnttS[dfRows[[1]][1]:dfRows[[1]][2],],
        currentEnttS[dfRows[[2]][1]:dfRows[[2]][2],]
    )
    dfFinal <- data.frame()
    colorList <- c("black", "red")
    for(i in seq(bc$tIDs)) {

        dfFinal <- rbind(dfFinal, data.frame(
            "X"=currenttS[[i]][,1],
            "Y"=currenttS[[i]][,2],
            "color"=colorList[i]
        ))
    }
    set.seed(4)
    shuffledDF <- dfFinal[sample(nrow(dfFinal)), ]
    enable("savetSNEPlot")
    plot(shuffledDF[[1]], shuffledDF[[2]], col=shuffledDF[[3]],
         xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=cex, lwd=0,
         xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
    axisTtSNE(xLim, yLim)
    leg <- rep(NA, length(listGroupORSamp))
    tempLegend <- legend("topright", bty="n",
                         inset=c(-0.25, -0.025), pch=15, lty=0, legend=leg,
                         cex=1.3, pt.cex=4, y.intersp=1.6,
                         text.width=strwidth("100"), xpd=TRUE,
                         lwd=2, pt.bg=colorList, col=colorList)
    multiplier <- 1.2
    strW <- strwidth(listGroupORSamp)*multiplier
    tX <- strW + tempLegend$text$x*1.01
    tY <- tempLegend$text$y
    text(tX, tY, listGroupORSamp, pos=2, xpd=NA, font=2, cex=multiplier)
}

tSGenHeat <- function(cex, xLim, yLim, tSHighl, dfRows, currenttS) {
    col1 <- colorRampPalette(c("#5E4FA2", "#3288BD"))(10)
    col2 <- colorRampPalette(c("#3288BD", "#ABDDA4"))(10)
    col3 <- rep("#ABDDA4", 5)
    col4 <- colorRampPalette(c("#ABDDA4", "#E6F598"))(10)
    col5 <- colorRampPalette(c("#E6F598", "#D53E4F"))(10)
    col <- c(col1, col2, col3, col4, col5)
    dataF <- list()
    for(i in seq(ncol(bc$concat))) {
        dataF[[i]] <- bc$concat[,i]
    }
    dataF <- data.frame(unlist(dataF))
    dataF <- cbind2(dataF, seq_len(nrow(dataF)))
    colnames(dataF) <- c("value", "ID")
    dataF <- dataF %>% arrange(dataF$value)
    dataF <- cbind2(dataF, colorRampPalette(col)(nrow(dataF)))
    dataF <- dataF %>% arrange(dataF$ID)
    colnames(dataF)[3] <- "col"
    legRange <- c(min(dataF[,1]), max(dataF[,1]))
    if(tSHighl == "") {
        highlight <- bc$availableparameters[[1]]
    } else {
        highlight <- tSHighl
    }
    highID <- which(colnames(bc$concat) == highlight)
    intStart <- (nrow(bc$concat)*(highID - 1)) + 1
    intEnd <- nrow(bc$concat)*highID
    dataF <- dataF[intStart:intEnd,]
    dataF["ID"] <- rownames(dataF) <- seq_len(nrow(dataF))
    enable("savetSNEPlot")
    cols <- dataF$col[dfRows[1]:dfRows[2]]
    plot(currenttS, xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=cex,
         lwd=0, col=cols, xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
    axisTtSNE(xLim, yLim)
    legendCol <- colorRampPalette(col)(80)
    len <- length(legendCol)
    barXCoords <- seq(xLim[1], xLim[2], length.out=len)
    y <- yLim[2]*1.15
    for(i in seq(legendCol)) {
        x <- barXCoords[i]*0.5
        points(x, y, pch=15, col=legendCol[i], xpd=TRUE, cex=2)
    }
    tX <- barXCoords[1]*0.6
    labs <- format(round(legRange[1], 2), nsmall=2)
    text(tX, y, xpd=TRUE, labels=labs, adj=1, cex=1.1)
    tX <- barXCoords[length(barXCoords)]*0.6
    labs <- format(round(legRange[2], 2), nsmall=2)
    text(tX, y, xpd=TRUE, labels=labs, adj=0, cex=1.1)
}

tSGenOvPo <- function(dfRows, currentEnttS, cex, xLim, yLim, tPops) {
    groupSampIndices <- seq(dfRows[1], dfRows[2])
    currenttS <- list()
    for(i in seq(tPops)) {
        id <- which(bc$concat[,tPops[i]] == 1)
        index <- intersect(groupSampIndices, id)
        currenttS[[i]] <- currentEnttS[index,]
    }
    dfFinal <- data.frame()
    colorList <- c("red2", "green4", "blue3", "orange", "black",
                   "wheat3", "brown", "orchid3", "salmon")
    for(i in seq(tPops)) {
        dfFinal <- rbind(dfFinal, data.frame("X"=currenttS[[i]][,1],
                                             "Y"=currenttS[[i]][,2],
                                             "color"=colorList[i]))
    }
    set.seed(5)
    shuffledDF <- dfFinal[sample(nrow(dfFinal)), ]
    enable("savetSNEPlot")
    plot(shuffledDF[[1]], shuffledDF[[2]], col=shuffledDF[[3]],
         xaxt="n", yaxt="n", ann=FALSE, pch=20, cex=cex, lwd=0,
         xlim=xLim, ylim=yLim, xaxs="i", yaxs="i")
    axisTtSNE(xLim, yLim)
    cols <- colorList[seq_len(length(tPops))]
    leg <- rep(NA, length(tPops))
    tempLegend <- legend("topright", bty="n", inset=c(-0.25, -0.025), pch=15,
                         lty=0, lwd=2, cex=1.3, pt.cex=4, xpd=TRUE,
                         y.intersp=1.6, text.width=strwidth("100"),
                         pt.bg=cols, col=cols, legend=leg)
    multiplier <- 1.2
    strW <- strwidth(tPops)
    tX <- strW*multiplier + tempLegend$text$x*1.01
    tY <- tempLegend$text$y
    text(tX, tY, tPops, pos=2, xpd=NA, font=2,
         cex=multiplier)
    hide("showingtSNEEvents")
}

showtSNE <- function() {
    hide("tSNEGroups")
    hide("tSNESamples")
    hide("tSPar")
    hide("tSNEParameters")
    hide("tSEvs")
    shinyjs::show("tSNEMode")
    shinyjs::show("tSNEDotSize")
    shinyjs::show("tSHighl")
    shinyjs::show("tSGroupOrSamp")
    shinyjs::show("tSGroupOrSampID")
    shinyjs::show("showingtSNEEvents")
    shinyjs::show("savetSNEPlot")
    enable("tSNEGenerate")
}

hidetSNE <- function() {
    hide("tSNEMode")
    hide("tSNEDotSize")
    hide("tSHighl")
    hide("tSGroupOrSamp")
    hide("tSGroupOrSampID")
    hide("tSGroupOrSampIDs")
    hide("tSNEPopulations")
    hide("savetSNEPlot")
    hide("showingtSNEEvents")
    shinyjs::show("tSNEGroups")
    shinyjs::show("tSNESamples")
    shinyjs::show("tSPar")
    shinyjs::show("tSNEParameters")
    shinyjs::show("tSEvs")
}

#UI----
ui <- fluidPage(
    useShinyjs(),
    extendShinyjs(text=jsCode, functions=c('disableTab', 'enableTab')),
    inlineCSS(css),
    tags$style("#about {border-color:white; font-size:12px}"),
    headerPanel(
        fluidRow(
            column(
                width=3,
                div(img(src=dataURI(
                    file=system.file("icons", "Logo.png", package="BCyto"),
                    mime="image/png"), height="50px"), "BCyto")
            ),
            column(
                width=2, offset=7, align="right",
                actionButton(inputId="about", label="About BCyto")
            )
        ),
    ),
    tags$style("#tabs {border-color:black; top:25px}"),
    tags$style(
        HTML(".tabbable > .nav > li > a {color:black}")
    ),
    tags$style(
        HTML(".tabbable > .nav > li[class=active] > a
                   {color:black; border-top-color:black;
                   border-left-color:black; border-right-color:black}")
    ),
    tags$style("#exportimagePlot {border-color:black}"),
    tags$style("#plotHelp {border-color:white}"),
    tags$style("#parentHelp {border-color:white}"),
    tags$style("#nextSample .fa{font-size: 11px}"),
    tags$style("#previousSample .fa{font-size: 11px}"),
    tags$style("#nextProlifS .fa{font-size: 11px}"),
    tags$style("#prevProlifS .fa{font-size: 11px}"),
    tabsetPanel(
        id="tabs",
        ##File----
        tabPanel(
            div(icon("cloud-upload-alt"), "File"),
            value="fileTab",
            br(),
            fluidRow(
                column(
                    width=4, offset=4,
                    sidebarPanel(
                        id="Fileleft", align="center",
                        shinyDirButton("directory", "Load FCS files",
                                       icon=icon("file-upload"),
                                       "Please select a folder with your
                                       FCS files"),
                        br(),
                        br(),
                        strong("or"),
                        br(),
                        br(),
                        shinyFilesButton("BCytoFileLoad", "Load BCyto file",
                                         icon=icon("file-upload"),
                                         "Please select a BCyto file (it
                                         needs to be located in the same
                                         folder of the FCS files)",
                                         multiple=FALSE),
                        br(),
                        br(),
                        strong("or"),
                        br(),
                        br(),
                        actionButton(inputId="loadTestData",
                                     label=div(icon("file-download"),
                                               "Download test data")),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        actionButton(inputId="saveFile", label=div(
                            icon("save"), "Save")),
                        br(),
                        br(),
                        tags$style("#fileHelp {border-color:white}"),
                        fluidRow(
                            column(
                                width=3, offset=9,
                                actionButton(inputId="fileHelp",
                                             icon("question-circle"))
                            )
                        ),
                        width=12
                    ),
                )
            ),
        ),
        ##Plot----
        tabPanel(
            div(
                icon("chart-area"), "Plot"), value="plotTab",
            br(),
            tags$style(".well {background-color:white; border-color:black}"),
            sidebarPanel(
                id="Plotleft",
                selectInput(inputId="samp", label=div(icon("vial"), "Sample"),
                            choices=c("")),
                tags$style("#customizeAxisY {border-color:black}"),
                fluidRow(
                    column(
                        width=8,
                        selectInput(inputId="Y",
                                    label=div(img(src=dataURI(
                                        file=system.file("icons", "Y.png",
                                                         package="BCyto"),
                                        mime="image/png"),
                                        height="15px"), "Y axis"),
                                    choices=c(""))
                    ),
                    column(
                        width=4, style="margin-top: 25px;",
                        actionButton(inputId="customizeAxisY", label="Scale")
                    )
                ),
                tags$style("#customizeAxisX {border-color:black}"),
                fluidRow(
                    column(
                        width=8,
                        selectInput(inputId="X",
                                    label=div(img(src=dataURI(
                                        file=system.file("icons", "X.png",
                                                         package="BCyto"),
                                        mime="image/png"),
                                        height="15px"), "X axis"),
                                    choices=c(""))
                    ),
                    column(
                        width=4, style="margin-top: 25px;",
                        actionButton(inputId="customizeAxisX", label="Scale")
                    )
                ),
                selectInput(inputId="typ", label= div(icon("chart-area"),
                                                      "Plot type"),
                            choices=c("Pseudocolor", "Contour", "Density",
                                      "Pseudocolor + Contour", "Histogram")),
                tags$style("#displayOptMain {border-color:black}"),
                fluidRow(
                    column(
                        width=4,
                        actionButton(inputId="displayOptMain",
                                     label=div(icon("bars"),
                                               "Display options"))
                    ),
                    column(
                        width=2, offset=5,
                        actionButton(inputId="plotHelp",
                                     icon("question-circle"))
                    )
                ),
                width=3
            ),
            mainPanel(
                id="Plotmiddle",
                fluidRow(
                    column(
                        width=3,
                        div(style="display:inline-block; position:fixed",
                            actionButton(inputId="previousSample",
                                         width="40px",
                                         icon("chevron-left"),
                                         style=bc$blackB),
                            actionButton(inputId="nextSample",
                                         width="40px",
                                         icon("chevron-right"),
                                         style=bc$blackB) )
                    ),
                    column(
                        width=5, offset=4, align="left", id="gatetools",
                        div(style="display:inline-block; position:fixed",
                            actionButton(inputId="rectang",
                                         width="40px",
                                         img(src=dataURI(
                                             file=system.file("icons",
                                                              "Rect.png",
                                                              package="BCyto"),
                                             mime="image/png"),
                                             height = "14px"),
                                         style=bc$blackB),
                            actionButton(inputId="polygon",
                                         width="40px",
                                         img(src=dataURI(
                                             file=system.file("icons",
                                                              "Polygon.png",
                                                              package="BCyto"),
                                             mime="image/png"),
                                             height = "14px"),
                                         style=bc$blackB),
                            actionButton(inputId="quad",
                                         width="40px",
                                         img(src=dataURI(
                                             file=system.file("icons",
                                                              "Quad.png",
                                                              package="BCyto"),
                                             mime="image/png"),
                                             height = "14px"),
                                         style=bc$blackB),
                            actionButton(inputId="interval",
                                         width="40px",
                                         img(src=dataURI(
                                             file=system.file("icons",
                                                              "Interval.png",
                                                              package="BCyto"),
                                             mime="image/png"),
                                             height = "7px"),
                                         style=bc$blackB),
                        )
                    ),
                    column(
                        width=5, offset=3,
                        div(style="display:inline-block; position:fixed",
                            textInput(inputId="drawGateName", label=NULL,
                                      placeholder="Please draw a gate",
                                      width="180px"), )
                    ),
                    column(
                        width=1, offset=0,
                        div(style="display:inline-block; position:fixed",
                            actionButton(inputId="gateOk", label="Ok",
                                         style=bc$blackB),
                            actionButton(inputId="gateCancel", label="Cancel",
                                         style=bc$whiteB)
                        )
                    ),
                ),
                br(),
                br(),
                uiOutput("plotUI"),
                tags$style("#saveMainPlot {border-color:black}"),
                div(style="display:inline-block",
                    downloadButton(outputId="saveMainPlot",
                                   label="Export image")),
                div(style="display:inline-block;vertical-align:top;
                    width: 20px;",
                    HTML("<br>")),
                div(style="display:inline-block", textOutput("events")),
                width=5
            ),
            sidebarPanel(
                id="Plotright",
                strong(div( icon("sitemap"), "Parent")),
                plotOutput("hierarchy", click="hierarchy_click",
                           dblclick="hierarchy_dblclick", height=450),
                fluidRow(
                    column(width=5,
                           tags$style("#exportImageGates {border-color:black}"
                           ),
                           downloadButton(outputId="exportImageGates",
                                          label="Export image"),
                    ),
                    column(width=4, offset=1,
                           tags$style("#editGate {border-color:black}"),
                           actionButton(inputId="editGate",
                                        label="Edit", icon("fas fa-pen")),
                    ),
                    column(width=1,
                           actionButton(inputId="parentHelp",
                                        icon("question-circle"))
                    )
                ),
                width=4
            )
        ),
        ##Compensation----
        tabPanel(
            div(icon("table"), "Compensation"), value="compTab",
            br(),
            sidebarPanel(
                selectInput(inputId="previewSample",
                            label=div(icon("vial"), "Preview sample"),
                            choices=c("")),
                checkboxInput(inputId="showUncomp",
                              label=strong("Overlay uncompensated"),
                              value=TRUE),
                selectInput(inputId="previewMatrix",
                            label=div(icon("table"), "Preview Matrix"),
                            choices=c("Cytometer-defined (applied)"=
                                          "Cytometer-defined")),
                tags$style("#createMatrix {border-color:black}"),
                actionButton(inputId="createMatrix",
                             label="Create New Matrix",
                             icon("plus-square")),
                br(),
                br(),
                tags$style(paste0("#saveMatrix {", bc$blackB, "}")),
                actionButton(inputId="saveMatrix",
                             label="Save New Matrix",
                             icon("sd-card")),
                br(),
                br(),
                tags$style("#cancelMatrix {border-color:black}"),
                actionButton(inputId="cancelMatrix",
                             label="Cancel"),
                br(),
                br(),
                tags$style("#applyMatrix {border-color:black}"),
                actionButton(inputId="applyMatrix",
                             label="Apply to Samples",
                             icon("share-square")),
                width=3
            ),
            mainPanel(
                fluidRow(
                    column(
                        width=12, rHandsontableOutput("comp")
                    )
                ),
                br(),
                br(),
                fluidRow(
                    column(
                        width=12, align="center", plotOutput("compGen")
                    )
                ),
                width=9
            )
        ),
        ##Ancestry----
        tabPanel(
            div(icon("chart-area"), "Ancestry plots"), value="ancestryTab",
            br(),
            sidebarPanel(
                selectInput(inputId="bgPreviewSample",
                            label=div(icon("vial"), "Sample"),
                            choices=c("")),
                selectInput(inputId="bgType",
                            label=div(icon("chart-area"), "Plot type"),
                            choices=c("", "Pseudocolor", "Contour", "Density",
                                      "Backgating"), selected=""),
                selectInput(inputId="bgPop",
                            label=div(tags$i(class="fas fa-clone"),
                                      HTML("Backgating <span style='color: red'
                                                      >Population</span>")),
                            choices=c("")),
                width=3
            ),
            mainPanel(
                plotOutput("ancestry", height=550),
                tags$style("#exportImageAncestry {border-color:black}"),
                downloadButton(outputId="exportImageAncestry",
                               label="Export image"),
                width=9
            )
        ),
        ##Overlays----
        tabPanel(
            div(icon("chart-area"), "Overlays"), value="overlayTab",
            br(),
            sidebarPanel(
                selectInput(inputId="ovTyp",
                            label=div(icon("chart-area"), "plot type"),
                            choices=c("", "Overlaid histogram",
                                      "Offset histogram", "Dot plot")),
                selectInput(inputId="ovTon",
                            label=div(icon("chart-area"), "Tone"),
                            choices=c("", "Colorful", "B&W")),
                selectInput(inputId="ovY",
                            label=div(img(src=dataURI(
                                file=system.file("icons", "Y.png",
                                                 package="BCyto"),
                                mime="image/png"), height="15px"),
                                "Y axis"), choices=c("")),
                selectInput(inputId="ovX",
                            label=div(img(src=dataURI(
                                file=system.file("icons", "X.png",
                                                 package="BCyto"),
                                mime="image/png"), height="15px"),
                                "X axis"), choices=c("")),
                selectInput(inputId="ovP",
                            label=div(icon("sitemap"), "Parent"),
                            choices=c("", "ungated")),
                tags$style("#ovSamples {border-color:black}"),
                actionButton(inputId="ovSamples",
                             label="Select samples", icon("plus-square")),
                br(),
                br(),
                tags$style("#displayOptOv {border-color:black}"),
                actionButton(inputId="displayOptOv",
                             label=div(icon("bars"), "Display options")),
                width=3
            ),
            mainPanel(
                plotOutput("overlays"),
                br(),
                br(),
                br(),
                tags$style("#ovSampleOrder {border-color:black}"),
                actionButton(inputId="ovSampleOrder",
                             label="Edit order", icon("fas fa-pen")),
                tags$style("#exportImageOverlay {border-color:black}"),
                downloadButton(outputId="exportImageOverlay",
                               label="Export image"),
                width=9
            ),
        ),
        ##Proliferation----
        tabPanel(
            div(icon("chart-area"), "Proliferation"), value="prolifTab",
            br(),
            sidebarPanel(
                tags$style("#step1 {border-color:black}"),
                actionButton(inputId="step1",
                             label=div("Step 1: gate undivided",
                                       tags$i(class="far fa-question-circle"))
                ),
                br(),
                br(),
                tags$style("#step2 {border-color:black}"),
                actionButton(inputId="step2",
                             label=div("Step 2: apply model",
                                       tags$i(class="far fa-question-circle"))
                ),
                br(),
                br(),
                tags$style("#step3 {border-color:black}"),
                actionButton(inputId="step3",
                             label=div("Step 3: results",
                                       tags$i(class="far fa-question-circle"))
                ),
                br(),
                br(),
                selectInput(inputId="prolifSButt",
                            label=div(icon("vial"), "Sample"),
                            choices=c("")),
                checkboxInput(inputId="prolifLabel", label="Numbers on top",
                              value=TRUE),
                checkboxInput(inputId="prolifGrid", label="Vertical lines",
                              value=TRUE),
                br(),
                br(),
                tags$style("#applyProlif {border-color:black}"),
                actionButton(inputId="applyProlif", label="Apply to Samples",
                             icon("share-square")),
                width=3
            ),
            mainPanel(
                htmlOutput("prolifStart"),
                fluidRow(
                    column(
                        width=3,
                        div(style="display:inline-block; position:fixed",
                            actionButton(inputId="prevProlifS", width="40px",
                                         icon("chevron-left"),
                                         style=bc$blackB),
                            actionButton(inputId="nextProlifS", width="40px",
                                         icon("chevron-right"),
                                         style=bc$blackB)
                        )
                    ),
                ),
                br(),
                br(),
                fluidRow(
                    column(
                        width=8,
                        plotOutput("prolifPlot", height=450),
                    ),
                    column(
                        width=4,
                        htmlOutput("prolifTable"),
                    ),
                ),
                fluidRow(
                    column(width=8,
                           tags$style("#exportImageProlif {border-color:black}"
                           ),
                           downloadButton(outputId="exportImageProlif",
                                          label="Export image"),
                    ),
                    column(width=4,
                           tags$style("#exportTableProlif {border-color:black}"
                           ),
                           downloadButton(outputId="exportTableProlif",
                                          label="Export complete table"),
                    ),
                ),
                width=9
            )
        ),
        ##t-SNE----
        tabPanel(
            div(icon("chart-pie"), "t-SNE"), value="tsneTab",
            br(),
            sidebarPanel(
                selectInput(inputId="tSNEMode",
                            label=div(icon("laptop-code"), "Mode"),
                            choices=c("Heatmap", "Overlay Groups or Samples",
                                      "Overlay Populations")),
                selectInput(inputId="tSHighl",
                            label=div(icon("highlighter"), "Highlight"),
                            choices=c("")),
                selectInput(inputId="tSGroupOrSamp",
                            label=div(icon("vial"), "Events to show"),
                            choices=c("All", "Group", "Sample")),
                selectInput(inputId="tSGroupOrSampID",
                            label=div(icon("vial"), "Group/Sample"),
                            choices=c("")),
                actionButton(inputId="tSGroupOrSampIDs",
                             label="Groups/Samples", icon("plus-square")),
                actionButton(inputId="tSNEPopulations",
                             label="Populations", icon("plus-square")),
                tags$style(HTML("[for=tSNEGroups]+span>.irs>.irs-single,
                        [for=tSNEGroups]+span>.irs-bar {background: black;
                        border-top: black; border-bottom: black}")),
                sliderInput("tSNEGroups", label="Number of groups",
                            min=1, max=6, step=1, value=1, ticks=FALSE),
                actionButton(inputId="tSNESamples", label="Select samples",
                             icon("plus-square")),
                br(),
                br(),
                tags$style(HTML("[for=tSNEDotSize]+span>.irs>.irs-single,
                        [for=tSNEDotSize]+span>.irs-bar {background: black;
                        border-top: black; border-bottom: black}")),
                sliderInput("tSNEDotSize", label="Dot size",
                            min=2, max=18, step=1, value=8, ticks=FALSE),
                selectInput(inputId="tSPar",
                            label=div(icon("sitemap"), "Parent"),
                            choices=c("", "ungated")),
                actionButton(inputId="tSNEParameters",
                             label="Select parameters", icon("plus-square")),
                br(),
                br(),
                selectInput(inputId="tSEvs",
                            label=div(icon("ellipsis-v"), "Events per sample"),
                            choices=c(" "=0, "1k"=1000, "5k"=5000, "10k"=10000,
                                      "25k"=25000, "50k"=50000)),
                br(),
                br(),
                actionButton(inputId="tSNEGenerate", label="Generate t-SNE",
                             icon("play-circle")),
                width=3
            ),
            mainPanel(
                plotOutput("tSNEPlot"),
                br(),
                br(),
                br(),
                tags$style("#savetSNEPlot {border-color:black}"),
                div(style="display:inline-block",
                    downloadButton(outputId="savetSNEPlot",
                                   label="Export image")),
                div(style="display:inline-block;vertical-align:top;
                    width: 20px;",
                    HTML("<br>")),
                div(style="display:inline-block",
                    textOutput("showingtSNEEvents")),
                width=9
            )
        ),
        ##Results----
        tabPanel(
            div(icon("table"), "Results"), value="resultTab",
            br(),
            sidebarPanel(
                selectInput(inputId="rsStat",
                            label=div(icon("percentage"), "Statistic"),
                            choices=c("", "Freq. of parent", "Freq. of...",
                                      "Freq. of total", "Median", "Count")),
                selectInput(inputId="rsParent",
                            label=div(icon("sitemap"),
                                      "Parent"), choices=c("")),
                selectInput(inputId="rsPop",
                            label=div(tags$i(class="fas fa-clone"),
                                      "Population"), choices=c("")),
                selectInput(inputId="rsCh",
                            label=div(
                                img(src=dataURI(
                                    file=system.file("icons", "X.png",
                                                     package="BCyto"),
                                    mime="image/png"),
                                    height="15px"),
                                "Parameter"),
                            choices=c("")),
                tags$style("#addResult {border-color:black}"),
                actionButton(inputId="addResult",
                             label="Add to table",
                             icon("plus-square")),
                tags$style("#editResult {border-color:black}"),
                width=3
            ),
            mainPanel(
                tags$head(tags$style(HTML("#results td, #results th
                                  {border-color:black}")) ),
                tableOutput("results"),
                actionButton(inputId="editResult",
                             label="Edit", icon("fas fa-pen")),
                tags$style("#exportTable {border-color:black}"),
                downloadButton(outputId="exportTable", label="Export table"),
                width=9
            )
        )
    )
)

#SERVER----
server <- function(input, output, session) {
    runjs("$('#drawGateName').attr('maxlength',18)")
    runjs("$('#drawGateName').bind('keypress', function(event) {
      var regex=new RegExp('^[/]+$');
      if(regex.test(String.fromCharCode(!event.charCode ? event.which:
      event.charCode))) {
          event.preventDefault()
      }
  })")
    hide("saveMatrix")
    hide("cancelMatrix")
    hide("prolifSButt")
    hide("prolifLabel")
    hide("prolifGrid")
    hide("prolifTable")
    hide("prevProlifS")
    hide("nextProlifS")
    hide("applyProlif")
    hide("exportImageProlif")
    hide("exportTableProlif")
    hidetSNE()
    disable("ovSamples")
    disable("step2")
    disable("step3")
    disable("exportTable")
    disable("editResult")
    disable("exportImageAncestry")
    disable("exportImageOverlay")
    disable("savetSNEPlot")
    disable("saveFile")
    js$disableTab("plotTab")
    js$disableTab("compTab")
    disableTabs()
    reactPar <- reactiveValues(d="root")
    reactNameChange <- reactiveValues(d=NULL)
    reactRectangle <- reactiveValues(d=FALSE)
    reactInterval <- reactiveValues(d=FALSE)
    reactPolygon <- reactiveValues(d=FALSE)
    reactQuadrant <- reactiveValues(d=FALSE)
    reactHover <- reactiveValues(d=0)
    hoverCoords <- reactiveValues(d=c(NULL))
    polygonCoords <- reactiveValues(d=data.frame(NULL))
    plotActivator <- reactiveValues(d=0)
    hierarchyActivator <- reactiveValues(d=0)
    bgActivator <- reactiveValues(d=0)
    reactGateFont <- reactiveValues(d=5)
    reactShowGateName <- reactiveValues(d=TRUE)
    reactAxisFont <- reactiveValues(d=15)
    reactShowAxis <- reactiveValues(d=TRUE)
    ovActivator <- reactiveValues(d=0)
    showOvPlot <- reactiveValues(d=FALSE)
    reactOvSamples <- reactiveValues(d=c())
    sampAlert <- reactiveValues(d="")
    reactModalTitle <- reactiveValues(d="")
    reactSampleOrder <- reactiveValues(d=NULL)
    compPlotActivator <- reactiveValues(d=0)
    reactReadOnly <- reactiveValues(d=TRUE)
    bc$appliedMatrix <- "Cytometer-defined"
    reactAxisCustom <- reactiveValues(d=NULL)
    prolifReady <- reactiveValues(d=FALSE)
    reactTSCh <- reactiveValues(d=NULL)
    reactTSSamp <- reactiveValues(d=NULL)
    reactGroupNames <- reactiveValues(d=NULL)
    reactTSPar <- reactiveValues(d=NULL)
    tSPlotActiv <- reactiveValues(d=0)
    reacttSIDs <- reactiveValues(d=NULL)
    reacttSNEPops <- reactiveValues(d=NULL)
    reactOvAxisFont <- reactiveValues(d=15)
    reactOvShowAxis <- reactiveValues(d=TRUE)
    reactDir <- reactiveValues(d=NULL)
    reactBCytoFile <- reactiveValues(d=NULL)
    reactHeight <- reactiveValues(d=200)
    reactOverw <- reactiveValues(d=NULL)
    reactTestData <- reactiveValues(d=FALSE)
    bc$currentTable <- NULL
    bc$entiretSNE <- NULL
    bc$betweenPeaks <- NULL
    bc$refPeaks <- NULL
    bc$adjustedL <- NULL
    bc$concatSamples <- NULL
    bc$tSNEListofGroups <- NULL
    bc$tSNEListofSamples <- NULL
    bc$loadingFile <- bc$ovBoolean <- FALSE

    #Plot----
    ##left----
    observeEvent(input$customizeAxisY, {
        reactAxisCustom$d <- input$Y
    })
    observeEvent(input$customizeAxisX, {
        reactAxisCustom$d <- input$X
    })

    observeEvent(reactAxisCustom$d, {
        reAddPops()
        index <- match(input$samp, bc$fileList)
        cs <- gs_pop_get_data(tempEnv$gs[index], reactPar$d)
        ff <- cytoframe_to_flowFrame(cs[[1]])@exprs[,reactAxisCustom$d]
        if(length(ff) > 0) {
            dataF <- data.frame(ff)
            if(nrow(dataF) > 5000) {
                set.seed(6)
                bc$tempDF <- dataF[sample.int(nrow(dataF), 5000),]
            } else {
                bc$tempDF <- dataF
            }
            value <- bc$customAxis[[which(bc$fluoCh == reactAxisCustom$d)]]
            showModal(modalDialog(
                tags$style(HTML(
                    paste0("[for=posslider]+span>.irs>.irs-single,
                                   [for=posslider]", bc$sliderCol))),
                tags$style(HTML(
                    paste0("[for=negslider]+span>.irs>.irs-single,
                                   [for=negslider]", bc$sliderCol))),
                tags$style(HTML(
                    paste0("[for=widthslider]+span>.irs>.irs-single,
                                   [for=widthslider]", bc$sliderCol))),
                fluidRow(
                    column(
                        width=9, align="right",
                        plotOutput("axisplot", height=440)),
                    column(
                        width=3, align="left", style="margin-top: 350px;",
                        actionButton(inputId="resetbutton", label="Reset",
                                     style=bc$whiteB))),
                fluidRow(
                    column(
                        width=12, align="center",
                        sliderInput("posslider", label="Positive Decades",
                                    min=2, max=7, step=0.02, value=value[1],
                                    ticks=FALSE),
                        sliderInput("negslider", label="Negative Decades",
                                    min=0, max=1, step=0.1, value=value[2],
                                    ticks=FALSE),
                        sliderInput("widthslider", label="Width Basis",
                                    min=0, max=3, step=0.1,
                                    value=log10(abs(value[3])),
                                    ticks=FALSE))),
                footer=list(actionButton(inputId="saveaxis", label="Ok",
                                         style=bc$blackB),
                            actionButton(inputId="cancelModal",
                                         label="Cancel", style=bc$whiteB)
                ),
                easyClose=FALSE,
                size="l"))
        } else {
            alert(paste0("The sample '", bc$fileList[index], "' has
            0 events in the selected parent. Please choose another sample or
                   change the parent."))
            reactAxisCustom$d <- NULL
        }
    })

    observeEvent(input$resetbutton, {
        updateSliderInput(inputId="posslider", value=4.42)
        updateSliderInput(inputId="negslider", value=0)
        updateSliderInput(inputId="widthslider", value=2)
    })

    output$axisplot <- renderPlot({
        par(mar=c(4,6,1,1) + 0.1, lwd=2)
        bc$widthBasis <- -10^input$widthslider
        if(input$widthslider < 1) {
            bc$widthBasis <- format(round(bc$widthBasis, 2), nsmall=2)
        } else if(input$widthslider < 2) {
            bc$widthBasis <- format(round(bc$widthBasis, 1), nsmall=1)
        } else if(input$widthslider <= 3) {
            bc$widthBasis <- round(bc$widthBasis)
        }
        trans <- flowjo_biexp(channelRange=4096, maxValue=262144,
                              pos=input$posslider,
                              neg=input$negslider,
                              widthBasis=as.numeric(bc$widthBasis))
        if(is(bc$tempDF, "data.frame")) {
            tempDF <- bc$tempDF[[1]]
        }
        tempDF <- trans(bc$tempDF)
        xLim <- c(0, 4100)
        if(length(tempDF) > 1) {
            referenceHist <- hist(tempDF, breaks=20, plot=FALSE)
            refCounts <- referenceHist$counts
            refDensity <- referenceHist$density
            multiplier <- refCounts/refDensity
            densLine <- density(tempDF)
            densLine$y <- densLine$y*multiplier[1]
            maxDens <- max(densLine[2]$y)
            plot(densLine, xaxt="n", yaxt="n", ann=FALSE, xlim=xLim,
                 ylim=c(0, maxDens*1.05), xaxs="i", yaxs="i")
            polygon(densLine, col="grey", border="black", lwd=1)
            axisTicks(2, "histo")
        } else {
            plot(0, xaxt="n", yaxt="n", ann=FALSE, xaxs="i", yaxs="i", cex=0)
        }
        title(ylab="Count", cex.lab=1+15/10, line=3, font.lab=2)
        index <- which(colnames(bc$gs) == reactAxisCustom$d)
        title(xlab=names(bc$ch)[index], cex.lab=1+15/10, line=3, font.lab=2)
        customAx <- c(input$posslider, input$negslider,
                      as.numeric(bc$widthBasis))
        axisTicks(1, "log", customAx)
    }, height=440, width=470)

    observeEvent(input$saveaxis, {
        currentA <- bc$customAxis[[which(bc$fluoCh == reactAxisCustom$d)]]
        newA <- c(input$posslider, input$negslider,
                  as.numeric(bc$widthBasis))
        if(!identical(currentA, newA)) {
            inverseT <- flowjo_biexp(channelRange=4096, maxValue=262144,
                                     pos=currentA[1],
                                     neg=currentA[2],
                                     widthBasis=currentA[3],
                                     inverse=TRUE)
            trans <- flowjo_biexp(channelRange=4096, maxValue=262144,
                                  pos=input$posslider,
                                  neg=input$negslider,
                                  widthBasis=as.numeric(bc$widthBasis))
            bc$customAxis[[which(bc$fluoCh == reactAxisCustom$d)]] <- newA
            popPaths <- gs_get_pop_paths(bc$gs)
            reAddPops()
            bc$gs <- gs_clone(tempEnv$gs)
            names(bc$customAxis) <- bc$fluoCh
            for(i in bc$fluoCh[bc$fluoCh != reactAxisCustom$d]) {
                transSub <- flowjo_biexp_trans(channelRange=4096,
                                               maxValue=262144,
                                               pos=bc$customAxis[[i]][1],
                                               neg=bc$customAxis[[i]][2],
                                               widthBasis=bc$customAxis[[i]][3]
                )
                bc$gs <- transform(bc$gs, transformerList(i, transSub))
            }
            transSub <- flowjo_biexp_trans(channelRange=4096,
                                           maxValue=262144,
                                           pos=input$posslider,
                                           neg=input$negslider,
                                           widthBasis=as.numeric(bc$widthBasis)
            )
            bc$gs <- transform(bc$gs, transformerList(reactAxisCustom$d,
                                                      transSub))
            recompute(bc$gs)
            for(i in seq(popPaths)[-1]) {
                gate <- gs_pop_get_gate(bc$gs, popPaths[i])[1][[1]]
                gateCh <- names(gate@parameters)
                for(j in seq(gateCh)) {
                    if(names(gate@parameters)[j] == reactAxisCustom$d) {
                        if(is(gate, "polygonGate")) {
                            toInvertObj <- gate@boundaries[,reactAxisCustom$d]
                            obj <- inverseT(toInvertObj)
                            gate@boundaries[,reactAxisCustom$d] <- trans(obj)
                        } else {
                            if(gate@min[[j]] != "-Inf") {
                                gate@min[[j]] <- trans(inverseT(gate@min[[j]]))
                            } else {
                                gate@min[[j]] <- -Inf
                            }
                            if(gate@max[[j]] != "Inf") {
                                gate@max[[j]] <- trans(inverseT(gate@max[[j]]))
                            } else {
                                gate@max[[j]] <- Inf
                            }
                        }
                        popgateList <- vapply(popgateList, function(x)
                            list(gate), list(length(popgateList)))
                        gs_pop_set_gate(bc$gs, gate@filterId, popgateList)
                        recompute(bc$gs, gate@filterId)
                    }
                }
            }
            reactAxisCustom$d <- NULL
            plotActivator$d <- plotActivator$d + 1
        } else {
            reactAxisCustom$d <- NULL
        }
        compPlotActivator$d <- compPlotActivator$d + 1
        hierarchyActivator$d <- hierarchyActivator$d + 1
        removeModal()
    })

    observeEvent(input$displayOptMain, {
        showModal(modalDialog(
            tags$style(HTML(paste0("[for=displayAxisFont]
            +span>.irs>.irs-single, [for=displayAxisFont]", bc$sliderCol))),
            checkboxInput("displayAxis", label="Show axis titles",
                          value=reactShowAxis$d),
            sliderInput("displayAxisFont", label=NULL, min=10, max=25,
                        value=reactAxisFont$d, ticks=FALSE),
            br(),
            tags$style(HTML(paste0("[for=displayGateFont]
            +span>.irs>.irs-single, [for=displayGateFont]", bc$sliderCol))),
            checkboxInput("displayGateName", label="Show gate name",
                          value=reactShowGateName$d),
            sliderInput("displayGateFont", label=NULL, min=1, max=30,
                        value=reactGateFont$d, ticks=FALSE),
            footer=list(
                actionButton(inputId="reloadplotmodal", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="s"))
        if(length(gs_get_pop_paths(bc$gs)) > 1) {
            shinyjs::show("displayGateName")
            shinyjs::show("displayGateFont")
        } else {
            hide("displayGateName")
            hide("displayGateFont")
        }
    })

    observeEvent(input$reloadplotmodal, {
        plotActivator$d <- plotActivator$d + 1
        removeModal()
        reactShowAxis$d <- input$displayAxis
        reactAxisFont$d <- input$displayAxisFont
        reactShowGateName$d <- input$displayGateName
        reactGateFont$d <- input$displayGateFont
    })

    observeEvent(c(input$displayAxis, input$displayOptMain) , {
        if(length(input$displayAxis) != 0) {
            if(input$displayAxis == TRUE) {
                shinyjs::show("displayAxisFont")
            } else {
                hide("displayAxisFont")
            }
        }
    })

    observeEvent(input$plotHelp, {
        showModal(modalDialog(
            title="Help",
            "Use the buttons at the top right side of the plot to
            interactively draw gates.",
            br(),
            br(),
            "Double clicking the population inside the gate leads to
            a new plot with the sorted population.",
            footer=NULL,
            size="w",
            easyClose=TRUE))
    })

    observeEvent(input$about, {
        showModal(modalDialog(
            title="BCyto (version 1.0)",
            "BCyto is an open-source project that provides an user-friendly,
            high-performance UI for Flow Cytometry analysis in R.",
            br(),
            br(),
            strong("Project website:"),
            "https://github.com/BonilhaCaio/BCyto",
            footer=NULL,
            size="m",
            easyClose=TRUE))
    })

    ##main----
    observeEvent(input$typ, {
        if(input$typ != "Histogram") {
            if(input$Y == "") {
                updateSelectInput(inputId="Y", selected=bc$Y)
            } else {
                plotActivator$d <- plotActivator$d + 1
            }
            enable("Y")
            disable("interval")
            enable("rectang")
            enable("polygon")
            enable("quad")
        } else {
            updateSelectInput(inputId="Y", selected="")
            disable("Y")
            enable("interval")
            disable("rectang")
            disable("polygon")
            disable("quad")
        }
    })

    observeEvent(input$nextSample, {
        current <- match(input$samp, bc$sampleList)
        if(current + 1 <= length(bc$sampleList)) {
            select <- bc$sampleList[[current + 1]]
            updateSelectInput(inputId="samp", selected=select)
        }
    })
    observeEvent(input$previousSample, {
        current <- match(input$samp, bc$sampleList)
        if(current > 1) {
            select <- bc$sampleList[[current - 1]]
            updateSelectInput(inputId="samp", selected=select)
        }
    })

    observeEvent(c(input$samp, input$nextSample, input$previousSample), {
        current <- match(input$samp, bc$sampleList)
        if(input$samp != "") {
            if(current >= length(bc$sampleList)) {
                disable("nextSample")
            } else {
                enable("nextSample")
            }
            if(current <= 1) {
                disable("previousSample")
            } else {
                enable("previousSample")
            }
        }
    })

    observeEvent(c(input$samp, input$Y, input$X, reactPar$d, input$typ,
                   input$gateCancel), {
                       reactRectangle$d <- FALSE
                       reactInterval$d <- FALSE
                       reactPolygon$d <- FALSE
                       reactQuadrant$d <- FALSE
                   }
    )

    observeEvent(input$gateOk, {
        "%notin%" <- Negate("%in%")
        if(input$drawGateName %notin% gs_get_pop_paths(bc$gs, path=1)) {
            if(reactRectangle$d == TRUE) {
                components <- c(input$main_brush$xmin,
                                input$main_brush$xmax,
                                input$main_brush$ymin,
                                input$main_brush$ymax)
                colnames <- c(input$X, input$Y)
                mat <- matrix(components, ncol=2,
                              dimnames=list(c("min", "max"), colnames))
                gate <- rectangleGate(.gate=mat)
            } else if(isolate(reactPolygon$d) == TRUE) {
                mat <- matrix(ncol=2, nrow=nrow(polygonCoords$d))
                for(i in seq(nrow(polygonCoords$d))) {
                    mat[i,] <- c(polygonCoords$d$x[i], polygonCoords$d$y[i])
                }
                colnames(mat) <- c(input$X, input$Y)
                gate <- polygonGate(mat)
            } else if(reactInterval$d == TRUE) {
                components <- c(input$main_brush$xmin,
                                input$main_brush$xmax)
                mat <- matrix(components, ncol=1,
                              dimnames=list(c("min", "max"), c(input$X)))
                gate <- rectangleGate(.gate=mat)
            }
            par <- isolate(reactPar$d)
            gs_pop_add(bc$gs, gate, parent=par, name=input$drawGateName)
            recompute(bc$gs)
            updateTextInput(inputId="drawGateName", value="")
            reactRectangle$d <- FALSE
            reactInterval$d <- FALSE
            reactPolygon$d <- FALSE
        } else {
            alert("Please choose a different gate name.")
        }
    })

    output$plotUI <- renderUI({
        if(input$typ == "Histogram") {
            if(reactInterval$d == TRUE) {
                brushOpts <- brushOpts("main_brush", fill="lightgrey",
                                       stroke="black", opacity=0.4, delay=10,
                                       direction="x")
            } else {
                brushOpts <- NULL
            }
            hoverOpts <- NULL
        } else {
            if(reactRectangle$d == TRUE) {
                brushOpts <- brushOpts("main_brush", fill="lightgrey",
                                       stroke="black", opacity=0.4, delay=10,
                                       direction="xy")
                hoverOpts <- NULL
            } else {
                brushOpts <- NULL
                if(isolate(reactPolygon$d) == TRUE
                   || reactQuadrant$d == TRUE) {
                    hoverOpts <- hoverOpts(id="mainplot_hover", delay=100,
                                           delayType="debounce",
                                           nullOutside=TRUE)
                } else {
                    hoverOpts <- NULL
                }
            }
        }
        plotOutput("mainplot", height=440, click="main_click",
                   dblclick="main_dblclick", brush=brushOpts, hover=hoverOpts)
    })

    output$mainplot <- renderPlot({
        par(mar=c(4,6,1,1) + 0.1, lwd=2)
        polygonCoords$d
        input$gateOk
        input$gateCancel
        plotActivator$d
        updateSelectInput(inputId="editGate", label="Edit")
        disable("editGate")
        reactNameChange$d <- NULL
        session$resetBrush("main_brush")
        ID <- match(input$samp, bc$fileList)
        par <- isolate(reactPar$d)
        newPlot(ID, input$X, input$Y, par, isolate(input$typ),
                isolate(reactShowAxis$d), isolate(reactAxisFont$d))
        detectGate(ID, input$X, input$Y, par, isolate(input$typ),
                   "plot", isolate(reactShowGateName$d),
                   isolate(reactGateFont$d))
        if(isolate(reactPolygon$d) == TRUE) {
            click <- isolate({input$main_click})
            if(!is.null(click)) {
                if(abs(click$x) > 0 && abs(click$y) > 0) {
                    coords <- polygonCoords$d
                    for(i in seq(nrow(coords))) {
                        if(i > 1) {
                            segments(coords$x[i-1], coords$y[i-1],
                                     coords$x[i], coords$y[i],
                                     col="red", lwd=2)
                        }
                        points(coords$x[i], coords$y[i], pch=19, cex=2,
                               col="red")
                    }
                }
            }
        } else if(reactQuadrant$d == TRUE) {
            abline(v=hoverCoords$d[1], h=hoverCoords$d[2], lwd=2, col="red")
        } else if(isolate(reactRectangle$d) == FALSE
                  && isolate(reactInterval$d) == FALSE
                  && isolate(reactPolygon$d) == FALSE
                  && isolate(reactQuadrant$d) == FALSE) {
            hide("gateOk", TRUE, "fade")
            hide("gateCancel", TRUE, "fade")
            hide("drawGateName", TRUE, "fade")
            delay(500, shinyjs::show("rectang", TRUE, "fade"))
            delay(500, shinyjs::show("polygon", TRUE, "fade"))
            delay(500, shinyjs::show("quad", TRUE, "fade"))
            delay(500, shinyjs::show("interval", TRUE, "fade"))
            updateTextInput(inputId="drawGateName", value= "",
                            placeholder="Please draw a gate")
            disable("drawGateName")
            polygonCoords$d <- data.frame(NULL)
        }
        output$events <- renderText({
            counts <- gs_pop_get_stats(bc$gs[input$samp])[,3]
            shownEv <- counts[which(gs_get_pop_paths(bc$gs) == par)][[1]]
            totalEv <- counts[1][[1]]
            paste0("Plotted events: ",
                   format(shownEv, big.mark=","), "/",
                   format(totalEv, big.mark=","))
        })
        if(grepl("SC", input$X) || input$X == "Time") {
            disable("customizeAxisX")
        } else {
            enable("customizeAxisX")
        }
        if(grepl("SC", input$Y) || input$Y == "Time"
           || isolate(input$typ) == "Histogram") {
            disable("customizeAxisY")
        } else {
            enable("customizeAxisY")
        }
    }, width=470)

    observeEvent(c(input$rectang, input$polygon, input$quad, input$interval), {
        if(!is.null(bc$gs)) {
            detectGate(1, input$X, input$Y, reactPar$d,
                       input$typ, ">= 4", isolate(reactShowGateName$d),
                       isolate(reactGateFont$d))
            if(length(bc$found) >= 4) {
                showModal(modalDialog(
                    "Only 4 gates are allowed on a plot. Would you like
                    to remove the 4 existing gates to draw new ones?",
                    footer=list(
                        actionButton(inputId="deleteAllModal",
                                     label="Delete existing gates",
                                     style=bc$blackB),
                        actionButton(inputId="cancelModal", label="Cancel",
                                     style=bc$whiteB)
                    ),
                    easyClose=FALSE,
                    size="s"))
            } else {
                hideTools()
            }
        }
    })

    observeEvent(input$rectang, {
        detectGate(1, input$X, input$Y, reactPar$d,
                   input$typ, ">= 4", isolate(reactShowGateName$d),
                   isolate(reactGateFont$d))
        if(length(bc$found) < 4) {
            reactRectangle$d <- TRUE
            secStep()
        }
    })

    observeEvent(input$polygon, {
        detectGate(1, input$X, input$Y, reactPar$d,
                   input$typ, ">= 4", isolate(reactShowGateName$d),
                   isolate(reactGateFont$d))
        if(length(bc$found) < 4) {
            reactPolygon$d <- TRUE
            secStep()
        }
    })

    observeEvent(input$quad, {
        detectGate(1, input$X, input$Y, reactPar$d,
                   input$typ, ">= 4", isolate(reactShowGateName$d),
                   isolate(reactGateFont$d))
        if(length(bc$found) < 4) {
            reactQuadrant$d <- TRUE
            delay(500, shinyjs::show("gateCancel", TRUE, "fade"))
        }
    })

    observeEvent(input$interval, {
        reactInterval$d <- TRUE
        hideTools()
        secStep()
    })

    observeEvent(input$mainplot_hover, {
        hover <- input$mainplot_hover
        coords <- hoverCoords$d
        if(!is.null(coords[1])) {
            if(max(hover$x, coords[1])/min(hover$x, coords[1]) > 1.1
               || max(hover$y, coords[2])/min(hover$y, coords[2]) > 1.1) {
                reactHover$d <- reactHover$d + 1
                hoverCoords$d[1] <- hover$x
                hoverCoords$d[2] <- hover$y
            }
        } else {
            reactHover$d <- reactHover$d + 1
            hoverCoords$d[1] <- hover$x
            hoverCoords$d[2] <- hover$y
        }
    })

    observeEvent(input$main_click, {
        inputX <- input$main_click$x
        inputY <- input$main_click$y
        if(reactPolygon$d == TRUE) {
            coords <- polygonCoords$d
            if(!is.null(coords)) {
                if(nrow(coords) > 2) {
                    maxX <- max(inputX, coords$x[1])
                    minX <- min(inputX, coords$x[1])
                    maxY <- max(inputY, coords$y[1])
                    minY <- min(inputY, coords$y[1])
                    if(grepl("SC", input$X)) {
                        xVerifier <- maxX/minX < 2.5
                    } else {
                        xVerifier <- maxX/minX < 1.1
                    }
                    if(grepl("SC", input$Y)) {
                        yVerifier <- maxY/minY < 2.5
                    } else {
                        yVerifier <- maxY/minY < 1.1
                    }
                    if(xVerifier && yVerifier) {
                        polygonCoords$d <- rbind(coords, c(coords$x[1],
                                                           coords$y[1]))
                        updateTextInput(inputId="drawGateName",
                                        placeholder="Type gate name")
                        enable("drawGateName")
                    } else {
                        nrowX <- coords$x[nrow(coords)]
                        nrowY <- coords$y[nrow(coords)]
                        if(nrowX != coords$x[1] && inputY != nrowY) {
                            polygonCoords$d <- rbind(coords, c(inputX, inputY))
                        }
                    }
                } else {
                    polygonCoords$d <- rbind(coords, c(inputX, inputY ))
                    colnames(polygonCoords$d) <- c("x", "y")
                }
            }
        } else if(reactQuadrant$d == TRUE) {
            mat <- matrix(c(inputX, inputY), ncol=2,
                          dimnames=list(c("value"),
                                        c(input$X, input$Y)))
            gate <- quadGate(.gate=mat)
            names <- c(paste0("Q1: ", input$X, "- ", input$Y, "+"),
                       paste0("Q2: ", input$X, "+ ", input$Y, "+"),
                       paste0("Q3: ", input$X, "+ ", input$Y, "-"),
                       paste0("Q4: ", input$X, "- ", input$Y, "-"))
            gs_pop_add(bc$gs, gate, parent=reactPar$d, names=names)
            recompute(bc$gs)
            detectGate(1, input$X, input$Y, reactPar$d,
                       input$typ, "OK", isolate(reactShowGateName$d),
                       isolate(reactGateFont$d))
            updateTextInput(inputId="drawGateName", value="")
            hierarchyActivator$d <- isolate(hierarchyActivator$d) + 1
            reactQuadrant$d <- FALSE
        }
    })

    observeEvent(c(input$main_brush, input$main_click, input$drawGateName), {
        if(reactRectangle$d == TRUE || reactInterval$d == TRUE) {
            if(is.null(input$main_brush)) {
                disable("gateOk")
                updateTextInput(inputId="drawGateName", value="",
                                placeholder="Please draw a gate")
                disable("drawGateName")
            } else {
                if(input$drawGateName == "") {
                    disable("gateOk")
                    updateTextInput(inputId="drawGateName",
                                    placeholder="Type gate name")
                    enable("drawGateName")
                } else {
                    enable("gateOk")
                }
            }
        }
        if(reactPolygon$d == TRUE) {
            if(input$drawGateName == "") {
                disable("gateOk")
            } else {
                enable("gateOk")
            }
        }
    })

    observeEvent(input$main_dblclick, {
        popPaths <- gs_get_pop_paths(bc$gs)
        if(length(popPaths) > 1) {
            inputX <- input$main_dblclick$x
            inputY <- input$main_dblclick$y
            popGates <- vapply(popPaths[-1], function(x)
                list(gs_pop_get_gate(bc$gs, x)[1][[1]]),
                list(length(popPaths[-1])))
            clickableGate <- list()
            for(i in seq(popGates)) {
                channels <- names(popGates[[i]]@parameters)
                if(length(channels) == 1) {
                    channels[2] <- bc$Y
                }
                popParent <- gs_pop_get_parent(bc$gs, popGates[[i]]@filterId)
                if(input$typ != "Histogram") {
                    if(channels[1] == input$X
                       && channels[2] == input$Y
                       && popParent == reactPar$d) {
                        clickableGate[[i]] <- popGates[[i]]
                    }
                } else {
                    if(channels[1] == input$X
                       && popParent == reactPar$d) {
                        clickableGate[[i]] <- popGates[[i]]
                    }
                }
            }
            updateParent <- FALSE
            singlePath <- gs_get_pop_paths(bc$gs, path=1)
            clickableGate <- unlist(clickableGate)
            for(i in clickableGate) {
                index <- which(singlePath == i@filterId)
                if(is(i, "rectangleGate")) {
                    if(length(names(i@parameters)) > 1) {
                        if(inputX >= i@min[[1]]
                           && inputX <= i@max[[1]]
                           && inputY >= i@min[[2]]
                           && inputY <= i@max[[2]]) {
                            updateParent <- TRUE
                            reactPar$d <- popPaths[index]
                        }
                    } else {
                        if(inputX >= i@min[[1]]
                           && inputX <= i@max[[1]]) {
                            updateParent <- TRUE
                            reactPar$d <- popPaths[index]
                        }
                    }
                }
                if(is(i, "polygonGate")) {
                    if(inputX >= min(i@boundaries[,1])
                       && inputX <= max(i@boundaries[,1])
                       && inputY >= min(i@boundaries[,2])
                       && inputY <= max(i@boundaries[,2])) {
                        updateParent <- TRUE
                        reactPar$d <- popPaths[index]
                    }
                }
            }
            if(updateParent == TRUE) {
                if(reactPar$d != "root") {
                    gate <- gs_pop_get_gate(bc$gs, reactPar$d)
                    channels <- names(gate[[1]]@parameters)
                }
                popParents <- gs_pop_get_count_fast(bc$gs[1], "freq")[,3][[1]]
                index <- which(popParents == reactPar$d)
                gatesWithThisParent <- popGates[index]
                if(length(gatesWithThisParent) > 0) {
                    channels <- names(gatesWithThisParent[[1]]@parameters)
                }
                updateSelectInput(inputId="X", selected=channels[1])
                if(length(channels) == 2) {
                    if(input$typ == "Histogram") {
                        updateSelectInput(inputId="typ",
                                          selected="Pseudocolor")
                    }
                    updateSelectInput(inputId="Y", selected=channels[2])
                    if(input$X == channels[1]
                       && input$Y == channels[2]) {
                        plotActivator$d <- plotActivator$d + 1
                    }
                } else {
                    if(input$typ != "Histogram") {
                        updateSelectInput(inputId="Y", selected="")
                        updateSelectInput(inputId="typ", selected="Histogram")
                    }
                    if(input$X == channels) {
                        plotActivator$d <- plotActivator$d + 1
                    }
                }
            }
        }
    })

    output$saveMainPlot <- downloadHandler(
        function() {
            paste0(substr(input$samp, 1, nchar(input$samp) - 4), ".png") },
        function(file) {
            png(file, units="in", height=6, width=6.44, res=300)
            par(mar=c(4,6,1,1) + 0.1, lwd=2)
            ID <- match(input$samp, bc$fileList)
            newPlot(ID, input$X, input$Y, reactPar$d,
                    input$typ, isolate(reactShowAxis$d),
                    isolate(reactAxisFont$d))
            detectGate(ID, input$X, input$Y, reactPar$d,
                       input$typ, "plot", isolate(reactShowGateName$d),
                       isolate(reactGateFont$d))
            dev.off()
        }
    )

    ##right----
    output$hierarchy <- renderPlot({
        input$gateOk
        hierarchyActivator$d
        input$deleteAllModal
        if(length(gs_get_pop_paths(bc$gs)) > 1) {
            hPlot <- plot(bc$gs)
            labels <- gs_get_pop_paths(bc$gs, path=1)
            labels[1] <- "ungated"
            for(i in seq(labels)) {
                if(substr(labels[i], 1, 1) == "Q" && nchar(labels[i]) >= 13) {
                    labels[i] <- substr(labels[i], 1, 2)
                }
            }
            names(labels) <- nodes(hPlot)
            nodeAttrs <- list(label=labels)
            attrs <- list(node=list(fillcolor="white", shape="box", width=1,
                                    color="gray90", style="rounded"),
                          graph=list(rankdir="TB"))
            index <- which(gs_get_pop_paths(bc$gs) == reactPar$d)
            if(is.null(reactNameChange$d)) {
                colour <- "black"
                names(colour) <- nodes(hPlot)[index]
            } else {
                if(reactPar$d == reactNameChange$d) {
                    colour <- "red"
                    names(colour) <- nodes(hPlot)[index]
                } else {
                    colour <- c("black", "red")
                    i2 <- which(gs_get_pop_paths(bc$gs) == reactNameChange$d)
                    names(colour) <- c(nodes(hPlot)[index], nodes(hPlot)[i2])
                }
            }
            nodeAttrs$color <- colour
            bc$hPlot <- plot(hPlot, nodeAttrs=nodeAttrs, attrs=attrs)
            plot(hPlot, nodeAttrs=nodeAttrs, attrs=attrs)
            js$enableTab("ancestryTab")
            js$enableTab("overlayTab")
            js$enableTab("prolifTab")
            js$enableTab("tsneTab")
            js$enableTab("resultTab")
        } else {
            nodes <- c("ungated", "B")
            edgeL <- list(ungated="B", B="ungated")
            graphNEL <- new("graphNEL", nodes=nodes, edgemode="directed",
                            edgeL=edgeL)
            attrs <- list(node=list(fillcolor="white", shape="rectangle",
                                    width=1))
            edgeAttrs <- list(color=c("ungated~B"="white"))
            nodeAttrs <- list(color=c("B"="white"), fontcolor=c("B"="white"))
            bc$hPlot <- plot(agopen(graphNEL, "", attrs=attrs,
                                    edgeAttrs=edgeAttrs, nodeAttrs=nodeAttrs))
            disableTabs()
        }
        pops <- gs_get_pop_paths(bc$gs, path=1)[-1]
        if(bc$loadingFile == FALSE) {
            updateSelectInput(inputId="bgPop", choices=c("", pops[-1]))
            updateSelectInput(inputId="ovP", choices=c("", pops))
        }
        updateSelectInput(inputId="tSPar", choices=c("", pops))
        updateSelectInput(inputId="rsParent", choices=c("", pops))
        bc$loadingFile <- FALSE
    })

    observeEvent(input$hierarchy_dblclick, {
        popPaths <- gs_get_pop_paths(bc$gs)
        if(length(popPaths) > 1) {
            agNode <- bc$hPlot@AgNode
            xy <- c(seq_len(length(agNode)))
            dataF <- data.frame(x=xy, y=xy)
            for(i in seq(agNode)) {
                rownames(dataF)[i] <- popPaths[i]
                dataF$x[i] <- agNode[[i]]@center@x
                dataF$y[i] <- agNode[[i]]@center@y
                dataF$short[i] <- agNode[[i]]@txtLabel@labelText
            }
            selected <- nearPoints(dataF, input$hierarchy_dblclick,
                                   xvar="x", yvar="y", threshold=25,
                                   maxpoints=1, addDist=TRUE)
            if(length(rownames(selected)) > 0) {
                popGates <- vapply(popPaths[-1], function(x)
                    list(gs_pop_get_gate(bc$gs, x)[1][[1]]),
                    list(length(popPaths[-1])))
                popParents <- gs_pop_get_count_fast(bc$gs[1], "freq")[,3][[1]]
                checkingPop <- rownames(selected)
                if(checkingPop != reactPar$d) {
                    reactPar$d <- checkingPop
                    if(reactPar$d != "root") {
                        gate <- gs_pop_get_gate(bc$gs, reactPar$d)
                        channels <- names(gate[[1]]@parameters)
                    }
                    index <- which(popParents == reactPar$d)
                    gatesWithThisParent <- popGates[index]
                    if(length(gatesWithThisParent) > 0) {
                        channels <- names(gatesWithThisParent[[1]]@parameters)
                    }
                    updateSelectInput(inputId="X", selected=channels[1])
                    if(length(channels) == 2) {
                        if(input$typ == "Histogram") {
                            updateSelectInput(inputId="typ",
                                              selected="Pseudocolor")
                        }
                        updateSelectInput(inputId="Y", selected=channels[2])
                        if(input$X == channels[1]
                           && input$Y == channels[2]) {
                            plotActivator$d <- plotActivator$d + 1
                        }
                    } else {
                        if(input$typ != "Histogram") {
                            updateSelectInput(inputId="Y",
                                              selected="")
                            updateSelectInput(inputId="typ",
                                              selected="Histogram")
                        }
                        if(input$X == channels) {
                            plotActivator$d <- plotActivator$d + 1
                        }
                    }
                }
            }
        }
    })

    observeEvent(input$hierarchy_click, {
        popPaths <- gs_get_pop_paths(bc$gs)
        agNode <- bc$hPlot@AgNode
        xy <- c(seq_len(length(agNode)))
        dataF <- data.frame(x=xy, y=xy)
        for(i in seq(agNode)) {
            rownames(dataF)[i] <- popPaths[i]
            dataF$x[i] <- agNode[[i]]@center@x
            dataF$y[i] <- agNode[[i]]@center@y
            dataF$short[i] <- agNode[[i]]@txtLabel@labelText
        }
        selected <- nearPoints(dataF, input$hierarchy_click,
                               xvar="x", yvar="y", threshold=25,
                               maxpoints=1, addDist=TRUE)
        if(length(rownames(selected)) > 0
           && rownames(selected) != "ungated") {
            reactNameChange$d <- rownames(selected)
            updateSelectInput(inputId="editGate",
                              label=HTML("<span style=
                                         'color: red'>Edit</span>"))
            enable("editGate")
        } else {
            updateSelectInput(inputId="editGate", label="Edit")
            disable("editGate")
            reactNameChange$d <- NULL
        }
    })

    observeEvent(input$editGate, {
        index <- which(gs_get_pop_paths(bc$gs) == reactNameChange$d)
        path <- gs_get_pop_paths(bc$gs, path=1)[index]
        showModal(modalDialog(
            if(substr(path, 1, 1) != "Q" || substr(path, 3, 4) != ": ") {
                paste0("Please type a new name for gate '", path, "':")
            } else{
                "Quadrant gates can't have their names changed."
            },
            br(),
            br(),
            if(substr(path, 1, 1) != "Q" || substr(path, 3, 4) != ": ") {
                textInput(inputId="newNameTextBox", label=NULL,
                          width="180px")
            },
            footer=list(
                if(substr(path, 1, 1) != "Q" || substr(path, 3, 4) != ": ") {
                    actionButton(inputId="changeNameModal", label="Ok",
                                 style=bc$blackB)
                },
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB),
                actionButton(inputId="deleteGate", label="Delete gate",
                             style=bc$blackB)
            ),
            easyClose=TRUE,
            size="m"))
    })

    observeEvent(input$deleteGate, {
        rootPlotLoad <- FALSE
        popPaths <- gs_get_pop_paths(bc$gs, path=1)
        index <- which(gs_get_pop_paths(bc$gs) == reactNameChange$d)
        gatetoDelete <- popPaths[index]
        i2 <- which(gs_get_pop_paths(bc$gs) == reactPar$d)
        showingParent <- popPaths[i2]
        popsAfter <- gh_pop_get_descendants(bc$gs[[1]], showingParent, path=1)
        popsAfterIDs <- vapply(popsAfter, function(x)
            which(popPaths == x), integer(1))
        popsBefore <- popPaths[-popsAfterIDs]
        gs_pop_remove(bc$gs, gatetoDelete)
        if(gatetoDelete == showingParent || gatetoDelete %in% popsBefore) {
            reactPar$d <- "root"
        }
        plotActivator$d <- plotActivator$d + 1
        hierarchyActivator$d <- isolate(hierarchyActivator$d) + 1
        removeModal()
        updateSelectInput(inputId="editGate",
                          label=HTML("<span style=
                                     'color: black'>Edit</span>"))
    })

    observeEvent(c(input$newNameTextBox, input$editGate), {
        if(!is.null(input$newNameTextBox)) {
            if(input$newNameTextBox == "") {
                disable("changeNameModal")
            } else {
                enable("changeNameModal")
            }
        }
    })

    observeEvent(input$changeNameModal, {
        popPaths <- gs_get_pop_paths(bc$gs, path=1)
        if(input$newNameTextBox %in% popPaths) {
            alert("Please choose a different gate name.")
        } else {
            index <- which(gs_get_pop_paths(bc$gs) == reactNameChange$d)
            gs_pop_set_name(bc$gs, popPaths[index], input$newNameTextBox)
            removeModal()
            if(reactNameChange$d == reactPar$d) {
                reactPar$d <- paste0("/", input$newNameTextBox)
            }
            plotActivator$d <- plotActivator$d + 1
        }
    })

    observeEvent(input$deleteAllModal, {
        for(i in bc$found) {
            gs_pop_remove(bc$gs, i)
        }
        removeModal()
        plotActivator$d <- plotActivator$d + 1
    })

    observeEvent(input$cancelModal, {
        reactAxisCustom$d <- NULL
        removeModal()
    })

    output$exportImageGates <- downloadHandler(
        "Gate hierarchy.png",
        function(file) {
            png(file, units="in", height=6, width=6.44, res=300)
            par(mar=c(4,6,1,1) + 0.1, lwd=2)
            plot(bc$hPlot)
            dev.off()
        }
    )

    observeEvent(input$parentHelp, {
        showModal(modalDialog(
            title="Help",
            "This sidebar shows the gates hierarchy.",
            br(),
            br(),
            "The current parent is highlighted by a black line.",
            br(),
            br(),
            "Double clicking a parent box leads to a new plot with the
            sorted population.",
            br(),
            br(),
            "Clicking in a parent box selects it to allow for editing.",
            footer=NULL,
            size="m",
            easyClose=TRUE))
    })

    #Compensation----
    ##left----
    observeEvent(input$createMatrix, {
        showModal(modalDialog(
            "Select the method.",
            br(),
            br(),
            if(is.null(bc$compDFs["AutoSpill"][[1]])) {
                actionButton(inputId="createAutoSpill", label="AutoSpill*",
                             style=bc$blackB)
            },
            actionButton(inputId="createmanually", label="Manually",
                         style=bc$blackB),
            br(),
            br(),
            if(is.null(bc$compDFs["AutoSpill"][[1]])) {
                "*Automated compensation calculation and matrix generation.
                This might take a few minutes"
            },
            footer= actionButton(inputId="cancelModal", label="Cancel",
                                 style=bc$whiteB),
            easyClose=FALSE,
            size="m"))
    })

    observeEvent(input$createAutoSpill, {
        if(!is.null(bc$compControlIDs)){
            onlyControls = bc$fileList[bc$compControlIDs]
            verifier = vapply(onlyControls, function(x)
                length(grep(x, list.files(bc$filePath))) > 0, TRUE)
            verifier = all(verifier)
            if(length(bc$compControlIDs) >= length(bc$fluoCh) && verifier) {
                hide("createmanually")
                hide("cancelModal")
                disable("createAutoSpill")
                message = "Please wait..."
                withProgress(message=message, detail="", value=0, max=100, {
                    controlDir = getwd()
                    controlDefFile = paste0(getwd(), "/fcs_control.csv")
                    unst = c(grep("Unstained", bc$fileList))
                    onlyControls = onlyControls[-unst]
                    cs = load_cytoset_from_fcs(bc$fileList)[onlyControls]
                    ffs = flowSet_to_list(cytoset_to_flowSet(cs))
                    tubeName = vapply(ffs, function(x)
                        x@description$`TUBE NAME`, character(1))
                    orderGuess = vapply(tubeName, function(x)
                        paste0(strsplit(x, " Stained Control")[[1]], "-A"),
                        character(1))
                    verifier = lapply(orderGuess, function(x)
                        any(x == bc$fluoCh))
                    verifier = all(unlist(verifier))
                    if(verifier == FALSE){
                        orderGuess = bc$fluoCh[order(bc$fluoCh)]
                    }
                    naVal = rep(NA, length(orderGuess))
                    futureCsv = data.frame("filename" = onlyControls,
                                           "dye" = orderGuess,
                                           "antigen" = naVal,
                                           "wavelength" = naVal)
                    write.table(futureCsv, file=controlDefFile, sep=",",
                                row.names=FALSE)
                    asp = get.autospill.param()
                    flowControl = suppressWarnings(
                        read.flow.control(controlDir, controlDefFile, asp)
                    )
                    flowGate = gate.flow.data(flowControl, asp)
                    spill = get.marker.spillover(TRUE, flowGate,
                                                 flowControl, asp)
                    refined = refine.spillover(spill, NULL, flowGate,
                                               flowControl, asp)
                })
                rightOrder = vapply(flowControl$marker.original, function(x)
                    which(x == bc$fluoCh), integer(1))
                reordered = cbind2(refined[[1]], rightOrder)
                reordered = reordered[order(reordered[,ncol(reordered)]),]
                reordered = reordered[,-ncol(reordered)]
                reordered = rbind2(reordered, rightOrder)
                reordered = reordered[,order(reordered[nrow(reordered),])]
                reordered = reordered[-nrow(reordered),]
                colnames(reordered) = rownames(reordered) = bc$fluoCh
                bc$compDFs["AutoSpill"][[1]] = reordered
                removeModal()
                choices = names(bc$compDFs)
                names(choices) = choices
                allNames = c()
                for(i in seq(choices)) {
                    if(choices[i] == bc$appliedMatrix) {
                        allNames[i] = paste(bc$appliedMatrix, "(applied)")
                    } else {
                        allNames[i] = choices[i]
                    }
                }
                names(choices) = allNames
                updateSelectInput(inputId="previewMatrix", choices=choices,
                                  selected="AutoSpill")
            } else {
                alert("AutoSpill requires all compensation control
                      files to be loaded.")
            }
        } else {
            alert("AutoSpill requires all compensation control
                  files to be loaded.")
        }
    })

    observeEvent(input$createmanually, {
        js$disableTab("fileTab")
        js$disableTab("plotTab")
        disableTabs()
        disable("createMatrix")
        disable("applyMatrix")
        shinyjs::show("saveMatrix")
        shinyjs::show("cancelMatrix")
        reactReadOnly$d <- FALSE
        removeModal()
    })

    observeEvent(input$cancelMatrix, {
        js$enableTab("fileTab")
        js$enableTab("plotTab")
        if(length(gs_get_pop_paths(bc$gs)) > 1) {
            js$enableTab("ancestryTab")
            js$enableTab("overlayTab")
            js$enableTab("prolifTab")
            js$enableTab("tsneTab")
            js$enableTab("resultTab")
        }
        hide("saveMatrix")
        hide("cancelMatrix")
        enable("createMatrix")
        reactReadOnly$d <- TRUE
        compPlotActivator$d <- compPlotActivator$d + 1
    })

    observeEvent(input$saveMatrix, {
        js$enableTab("fileTab")
        js$enableTab("plotTab")
        if(length(gs_get_pop_paths(bc$gs)) > 1) {
            js$enableTab("ancestryTab")
            js$enableTab("overlayTab")
            js$enableTab("prolifTab")
            js$enableTab("tsneTab")
            js$enableTab("resultTab")
        }
        hide("saveMatrix")
        hide("cancelMatrix")
        enable("createMatrix")
        reactReadOnly$d <- TRUE
        index <- length(grep("Manual Matrix", names(bc$compDFs))) + 1
        matName <- paste("Manual Matrix", index)
        bc$compDFs[matName][[1]] <- bc$currentTable/100
        matChoices()
        updateSelectInput(inputId="previewMatrix", choices=bc$choices,
                          selected=names(bc$compDFs)[length(bc$compDFs)])
    })

    observeEvent(input$applyMatrix, {
        showModal(modalDialog(
            paste0("The matrix '", input$previewMatrix,
                   "' will be applied to all samples.
                   Do you want to proceed?"),
            footer=list(
                actionButton(inputId="OKapplymatrix", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="s"))
    })

    observeEvent(input$OKapplymatrix, {
        plotActivator$d <- plotActivator$d + 1
        bc$appliedMatrix <- input$previewMatrix
        popPaths <- gs_get_pop_paths(bc$gs)
        reAddPops()
        bc$gs <- gs_clone(tempEnv$gs)
        names(bc$customAxis) <- bc$fluoCh
        for(i in seq(bc$fluoCh)) {
            trans <- flowjo_biexp_trans(channelRange=4096, maxValue=262144,
                                        pos=bc$customAxis[[i]][1],
                                        neg=bc$customAxis[[i]][2],
                                        widthBasis=bc$customAxis[[i]][3])
            bc$gs <- transform(bc$gs, transformerList(bc$fluoCh[i], trans))
        }
        recompute(bc$gs)
        matChoices()
        updateSelectInput(inputId="previewMatrix", choices=bc$choices,
                          selected=bc$appliedMatrix)
        removeModal()
        plotActivator$d
    })

    observeEvent(c(input$previewMatrix, input$cancelMatrix), {
        if(bc$appliedMatrix == input$previewMatrix) {
            disable("applyMatrix")
        } else {
            enable("applyMatrix")
        }
    })

    ##main----
    observeEvent(input$comp, {
        handsonobject <- isolate(input$comp)
        invertedtable <- unlist(handsonobject$d)
        table <- vector(mode="list", length=length(bc$fluoCh))
        index <- 0
        for(i in seq(bc$fluoCh)) {
            for(j in seq(bc$fluoCh)) {
                index <- index + 1
                table[[j]][i] <- invertedtable[[index]]
            }
        }
        table <- as.data.frame(table, row.names=bc$fluoCh)
        colnames(table) <- rownames(table)
        if(!is.null(bc$currentTable)
           && !identical(bc$currentTable, table)
           && reactReadOnly$d == FALSE) {
            compPlotActivator$d <- compPlotActivator$d + 1
        }
        bc$currentTable <- table
    })

    observe({
        output$compGen <- renderPlot({
            input$previewMatrix
            compPlotActivator$d
            ID <- match(input$previewSample, bc$fileList)
            if(isolate(reactReadOnly$d) == TRUE) {
                withProgress(message="Please wait...", detail="", value=0,
                             max=100, {
                                 compGen(ID,
                                         bc$compDFs[input$previewMatrix][[1]],
                                         input$showUncomp)
                             })
            } else {
                withProgress(message="Please wait...", detail="", value=0,
                             max=100, {
                                 compGen(ID,
                                         bc$currentTable/100,
                                         input$showUncomp)
                             })
            }
        }, height=reactHeight$d, width=reactHeight$d)
    })

    #Ancestry----
    ##left----
    observeEvent(input$bgType, {
        if(isolate(input$bgType) == "Backgating"
           && length(gs_get_pop_paths(bc$gs)) > 2) {
            enable("bgPop")
        } else {
            updateSelectInput(inputId="bgPop", selected="")
            disable("bgPop")
        }
    })

    observeEvent(input$bgPop, {
        if(input$bgType == "Backgating") {
            bgActivator$d <- isolate(bgActivator$d) + 1
        }
    })

    ##main----
    output$ancestry <- renderPlot({
        bgActivator$d
        ancestryGen(input$bgType, input$bgPreviewSample, isolate(input$bgPop))
    }, height=480, width=800)

    output$exportImageAncestry <- downloadHandler(
        "Ancestry plots.png",
        function(file) {
            png(file, units="in", height=6.6, width=11.05, res=300)
            ancestryGen(input$bgType, input$bgPreviewSample, input$bgPop)
            dev.off()
        }
    )

    #Overlays----
    ##left----
    observeEvent(c(input$ovTyp, input$ovTon, input$ovY, input$ovX,
                   input$ovP), {
                       samples <- isolate(reactOvSamples$d)
                       if(input$ovTyp == "Overlaid histogram"
                          || input$ovTyp == "Offset histogram") {
                           if(input$ovY != "") {
                               updateSelectInput(inputId="ovY", selected="")
                           } else {
                               ovActivator$d <- ovActivator$d + 1
                           }
                           disable("ovY")
                           if(input$ovTyp == "Overlaid histogram") {
                               if(input$ovX != "" && input$ovP != "") {
                                   enable("ovSamples")
                                   if(length(samples) > 2) {
                                       reactOvSamples$d <- samples[seq_len(3)]
                                   }
                                   if(length(samples) < 3) {
                                       showOvPlot$d <- TRUE
                                   } else {
                                       showOvPlot$d <- FALSE
                                   }
                                   showOvPlot$d <- TRUE
                               } else {
                                   disable("ovSamples")
                                   showOvPlot$d <- FALSE }
                               if(input$ovTon != "") {
                                   updateSelectInput(inputId="ovTon",
                                                     selected="")
                               } else {
                                   ovActivator$d <- ovActivator$d + 1
                               }
                               disable("ovTon")
                           } else {
                               if(input$ovTon != ""
                                  && input$ovX != ""
                                  && input$ovP != "") {
                                   enable("ovSamples")
                                   ovActivator$d <- ovActivator$d + 1
                                   showOvPlot$d <- TRUE
                               } else {
                                   disable("ovSamples")
                                   showOvPlot$d <- FALSE }
                               enable("ovTon")
                               ovActivator$d <- ovActivator$d + 1
                           }
                       } else {
                           if(length(samples) > 2) {
                               reactOvSamples$d <- samples[seq_len(2)]
                           }
                           if(input$ovTyp != ""
                              && input$ovTon != ""
                              && input$ovY != ""
                              && input$ovX != ""
                              && input$ovP != "") {
                               enable("ovSamples")
                               ovActivator$d <- ovActivator$d + 1
                               showOvPlot$d <- TRUE
                           } else {
                               disable("ovSamples")
                               showOvPlot$d <- FALSE
                           }
                           enable("ovTon")
                           enable("ovY")
                           ovActivator$d <- ovActivator$d + 1
                       }
                   }
    )

    observeEvent(input$ovSamples, {
        sampAlert$d <- ""
        updateCheckboxGroupInput(inputId="samplecheckbox",
                                 selected=isolate(reactOvSamples$d))
        if(input$ovTyp == "Dot plot") {
            modaltitle <- "Select 2 samples for the Dot plot."
        } else if(input$ovTyp == "Overlaid histogram") {
            modaltitle <- "Select 2-3 samples for the Overlaid histogram."
        } else {
            modaltitle <- "Select 2-10 samples for the Offset histogram."
        }
        showModal(modalDialog(
            tags$style("#selectedsamplealert{color: red}"),
            strong(modaltitle),
            checkboxGroupInput("samplecheckbox", label="",
                               choices=bc$sampleList, width="100%" ),
            footer=fluidRow(
                column(width=8, align="left",
                       textOutput("selectedsamplealert")
                ),
                column(width=4, align="right", list(
                    actionButton(inputId="ovokbutton", label="Ok",
                                 style=bc$blackB),
                    actionButton(inputId="cancelModal", label="Cancel",
                                 style=bc$whiteB)
                ),
                )
            ),
            easyClose=TRUE,
            size="m"))
        disable("ovokbutton")
    })

    observeEvent(input$samplecheckbox, {
        len <- length(input$samplecheckbox)
        if(input$ovTyp == "Dot plot") {
            if(len == 2) {
                sampAlert$d <- ""
                enable("ovokbutton")
            } else {
                if(len > 2) {
                    sampAlert$d <- "Only 2 samples are allowed for Dot plots."
                } else {
                    sampAlert$d <- ""
                }
                disable("ovokbutton") }
        } else if(input$ovTyp == "Overlaid histogram") {
            sampAlert$d <- "Up to 3 samples are allowed for
            Overlaid histograms."
            if(len == 2 || len == 3) {
                sampAlert$d <- ""
                enable("ovokbutton")
            } else {
                if(len > 3) {
                    sampAlert$d <- "Up to 3 samples are allowed for
                    Overlaid histograms."
                } else {
                    sampAlert$d <- "" }
                disable("ovokbutton")
            }
        } else {
            sampAlert$d <- "Up to 10 samples are allowed for
            Overlaid histograms."
            if(len > 1 && len < 11) {
                sampAlert$d <- ""
                enable("ovokbutton")
            } else {
                if(len > 10) {
                    sampAlert$d <- "Up to 10 samples are allowed for
                    Overlaid histograms."
                } else {
                    sampAlert$d <- ""
                }
                disable("ovokbutton")
            }
        }
    })

    output$selectedsamplealert <- renderText({
        sampAlert$d
    })

    observeEvent(input$ovokbutton, {
        ovActivator$d <- ovActivator$d + 1
        reactOvSamples$d <- input$samplecheckbox
        removeModal()
    })

    observeEvent(input$displayOptOv, {
        t1 <- "[for=ovdisplayaxisfont]+span>.irs>.irs-single, "
        t2 <- "[for=ovdisplayaxisfont]"
        showModal(modalDialog(
            tags$style(HTML(paste0(t1, t2, bc$sliderCol))),
            checkboxInput("ovdisplayaxis", label="Show axis titles",
                          value=reactOvShowAxis$d),
            sliderInput("ovdisplayaxisfont", label=NULL, min=10, max=25,
                        value=reactOvAxisFont$d, ticks=FALSE),
            footer=list(
                actionButton(inputId="ovreloadplotmodal", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="s"))
    })

    observeEvent(input$ovreloadplotmodal, {
        ovActivator$d <- ovActivator$d + 1
        removeModal()
        reactOvShowAxis$d <- input$ovdisplayaxis
        reactOvAxisFont$d <- input$ovdisplayaxisfont
    })

    observeEvent(c(input$ovdisplayaxis, input$displayOptOv) , {
        if(length(input$ovdisplayaxis) != 0) {
            if(input$ovdisplayaxis == TRUE) {
                shinyjs::show("ovdisplayaxisfont")
            } else {
                hide("ovdisplayaxisfont")
            }
        }
    })

    ##main----
    output$overlays <- renderPlot({
        ovActivator$d
        if(!is.null(reactOvSamples$d)) {
            if(isolate(showOvPlot$d) == TRUE) {
                enable("displayOptOv")
                enable("ovSampleOrder")
                enable("exportImageOverlay")
                if(isolate(input$ovP) == "ungated") {
                    currentParent <- "root"
                } else {
                    currentParent <- isolate(input$ovP)
                }
                bc$ovBoolean <- TRUE
                ID <- match(isolate(reactOvSamples$d), bc$fileList)
                overlay(ID, isolate(input$ovX), isolate(input$ovY),
                        currentParent, isolate(input$ovTyp),
                        isolate(input$ovTon), isolate(reactOvShowAxis$d),
                        isolate(reactOvAxisFont$d))
            } else {
                bc$ovBoolean <- FALSE
                disable("displayOptOv")
                disable("ovSampleOrder")
                disable("exportImageOverlay")
            }
        } else {
            bc$ovBoolean <- FALSE
            disable("displayOptOv")
            disable("ovSampleOrder")
            disable("exportImageOverlay")
        }
    }, height=423, width=800)

    observeEvent(input$ovSampleOrder, {
        reactModalTitle$d <- 1
        reactSampleOrder$d <- NULL
        choices <-
            showModal(modalDialog(
                strong(textOutput("ordermodaltitle")),
                radioButtons("orderradiobutton", label="",
                             choices=rev(isolate(reactOvSamples$d)),
                             width="100%",
                             selected="" ),
                footer=actionButton(inputId="cancelModal", label="Cancel",
                                    style=bc$whiteB),
                easyClose=TRUE,
                size="m"))
    })

    output$ordermodaltitle <- renderText({
        paste0("Select sample number ", reactModalTitle$d,
               " (from top to bottom).")
    })

    observeEvent(input$orderradiobutton, {
        reactModalTitle$d <- reactModalTitle$d + 1
        ID <- match(input$orderradiobutton, reactOvSamples$d)
        for(i in seq(input$orderradiobutton)) {
            reactSampleOrder$d <- append(reactSampleOrder$d, ID)
        }
        if(length(reactSampleOrder$d) != length(reactOvSamples$d)) {
            samplestoShow <- reactOvSamples$d
            for(i in rev(sort(reactSampleOrder$d))) {
                preShow <- samplestoShow != samplestoShow[i]
                samplestoShow <- samplestoShow[preShow]
            }
            updateRadioButtons(inputId="orderradiobutton",
                               choices=rev(samplestoShow), selected="")
        } else {
            reactOvSamples$d <- rev(reactOvSamples$d[reactSampleOrder$d])
            ovActivator$d <- ovActivator$d + 1
            removeModal()
        }
    })

    output$exportImageOverlay <- downloadHandler(
        function() {
            paste0(input$ovTyp, ".png")
        },
        function(file) {
            png(file, units="in", height=5.8, width=11.05, res=300)
            if(isolate(input$ovP) == "ungated") {
                currentParent <- "root"
            } else {
                currentParent <- isolate(input$ovP)
            }
            ID <- match(isolate(reactOvSamples$d), bc$fileList)
            overlay(ID, isolate(input$ovX), isolate(input$ovY), currentParent,
                    isolate(input$ovTyp), isolate(input$ovTon),
                    isolate(reactOvShowAxis$d), isolate(reactOvAxisFont$d))
            dev.off()
        }
    )

    #Proliferation----
    ##left----
    observeEvent(input$applyProlif, {
        prolifReady$d <- TRUE
    })

    observeEvent(input$step1, {
        showModal(modalDialog(
            title = "Step 1 help",
            "Follow the instructions to start the proliferation tool.",
            footer=NULL,
            easyClose=TRUE,
            size="s"))
    })
    observeEvent(input$step2, {
        showModal(modalDialog(
            title = "Step 2 help",
            "Navigate through samples and choose the best fitting model
            to apply.",
            footer=NULL,
            easyClose=TRUE,
            size="s"))
    })
    observeEvent(input$step3, {
        showModal(modalDialog(
            title = "Step 3 help",
            "You can use the generated results or return to step 2
            to apply a new model.",
            footer=list(
                actionButton(inputId="goToStep2", label="Apply a new model",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="s"))
    })
    observeEvent(input$goToStep2, {
        prolifReady$d <- FALSE
        removeModal()
    })

    ##main----
    output$prolifStart <- renderUI({
        l0 <- "This tool automates the detection and quantification of
        division peaks from cell proliferation assays."
        l1 <- "To initialize the proliferation tool:"
        l2 <- "-go to the Plot tab;"
        l3 <- "-select the desired parent;"
        l4 <- "-select the proliferation dye in the X axis
        (i.e. FITC-A :: CFSE);"
        l5 <- "-select 'Histogram' in 'Plot type';"
        l6 <- "-make a narrow gate (that does not reach the X axis limit)
        on any sample for undivided cells (brightest peak);"
        l7 <- "-certify that the population is named as 'undivided';"
        l8 <- "-return to this tab."
        HTML("<br/>", paste(strong(l0)), "<br/>", "<br/>",
             paste(strong(l1),l2,l3,l4,l5,l6,l7,l8, sep="<br/>"))
    })

    observeEvent(input$nextProlifS, {
        ID <- match(input$prolifSButt, bc$sampleList)
        if(ID + 1 <= length(bc$sampleList)) {
            updateSelectInput(inputId="prolifSButt",
                              selected=bc$sampleList[[ID + 1]])
        }
    })
    observeEvent(input$prevProlifS, {
        ID <- match(input$prolifSButt, bc$sampleList)
        if(ID > 1) {
            updateSelectInput(inputId="prolifSButt",
                              selected=bc$sampleList[[ID - 1]])
        }
    })
    observeEvent(c(input$prolifSButt, input$nextProlifS, input$prevProlifS), {
        if(input$prolifSButt != "") {
            ID <- match(input$prolifSButt, bc$sampleList)
            if(ID >= length(bc$sampleList)) {
                disable("nextProlifS")
            } else {
                enable("nextProlifS")
            }
            if(ID <= 1) {
                disable("prevProlifS")
            } else {
                enable("prevProlifS")
            }
        }
    })

    output$prolifPlot <- renderPlot({
        input$tabs
        popPaths <- gs_get_pop_paths(bc$gs, path=1)
        if("undivided" %in% popPaths) {
            prolifGen(input$prolifSButt, isolate(reactShowAxis$d),
                      isolate(reactAxisFont$d), prolifReady$d,
                      input$prolifLabel, input$prolifGrid)
        } else {
            prolifOff()
            hide("prolifPlot")
            hide("applyProlif")
            disable("step2")
            hide("prolifSButt")
            hide("prevProlifS")
            hide("nextProlifS")
            shinyjs::show("prolifStart")
            enable("step1")
        }
    }, height=440, width=470)

    output$exportImageProlif <- downloadHandler(
        function() {
            paste0(substr(input$prolifSButt, 1, nchar(input$prolifSButt) - 4),
                   ".png")
        },
        function(file) {
            png(file, units="in", height=6, width=6.44, res=300)
            popPaths= gs_get_pop_paths(bc$gs, path=1)
            if("undivided" %in% popPaths) {
                prolifGen(input$prolifSButt, isolate(reactShowAxis$d),
                          isolate(reactAxisFont$d), prolifReady$d,
                          input$prolifLabel, input$prolifGrid)
            }
            dev.off()
        }
    )

    ##right----
    output$prolifTable <- renderUI({
        input$tabs
        #input$prolifSButt
        if(prolifReady$d == TRUE) {
            prolifTableGen(input$prolifSButt)
        }
    })

    output$exportTableProlif <- downloadHandler(
        "Proliferation results.csv",
        function(file) {
            colNames = c("Number of divisions", "Divided cells",
                         "Undivided cells")
            for(i in 2:nrow(bc$refPeaks)) {
                colNames[i+2] = paste0("Division ", (i - 1))
            }
            colNames[length(colNames) + 1] = "Percent divided"
            colNames[length(colNames) + 1] = "Percent undivided"
            constantColLenght = length(colNames)
            for(i in 2:nrow(bc$refPeaks)) {
                pasteString = paste0("Percent division ", (i - 1))
                colNames[i+(constantColLenght-1)] = pasteString
            }
            completeTable = data.frame(matrix("", ncol=length(colNames),
                                              nrow=length(bc$sampleList)))
            for(j in seq(bc$sampleList)) {
                cs = gs_pop_get_data(bc$gs[bc$sampleList[j]], bc$prolifParent)
                df = as.data.frame.array(cytoframe_to_flowFrame(cs[[1]])@exprs)
                dfX = df[,bc$prolifChannel]
                bc$divPercents = c()
                bc$divCounts = c()
                bc$total = length(dfX)
                for(i in seq(bc$gatecoords)) {
                    filter = dfX[dfX > bc$gatecoords[[i]][1]]
                    filter = filter[filter < bc$gatecoords[[i]][2]]
                    bc$divCounts[i] = length(filter)
                    bc$divPercents[i] = length(filter)*100/bc$total
                }
                bc$divCounts = rev(bc$divCounts)
                round = round(bc$divPercents, 2)
                bc$divPercents = rev(as.numeric(format(round, nsmall=2)))
                bc$totalDivCount = bc$total-bc$divCounts[1]
                bc$totalDivPerc = 100-bc$divPercents[1]
                columnValues = c((nrow(bc$refPeaks) - 1),
                                 bc$totalDivCount, bc$divCounts[1])
                for(i in 2:nrow(bc$refPeaks)) {
                    columnValues[i+2] = bc$divCounts[i]
                }
                ID = length(columnValues) + 1
                round = round(bc$totalDivPerc, 2)
                columnValues[ID] = as.numeric(format(round, nsmall=2))
                columnValues[ID] = bc$divPercents[1]
                constantColLenght = length(columnValues)
                for(i in 2:nrow(bc$refPeaks)) {
                    columnValues[i+(constantColLenght-1)] = bc$divPercents[i]
                }
                completeTable[j,] = columnValues
            }
            colnames(completeTable) = colNames
            rownames(completeTable) = bc$sampleList
            write.csv(completeTable, file)
        }
    )

    #t-SNE----
    ##left----
    observeEvent(input$tSNEGroups, {
        if(length(reactTSSamp$d) != input$tSNEGroups) {
            reactTSSamp$d <- NULL
        }
    })

    observeEvent(input$tSNESamples, {
        sampleCol <- list()
        for(i in seq(bc$sampleList)) {
            sampleCol[[i]] <- list(br(), br(), bc$sampleList[i])
        }
        groupCols <- list()
        for(i in seq(input$tSNEGroups)) {
            groupCols[[i]] <- list()
        }
        if(!is.null(reactTSSamp$d)) {
            for(j in seq(input$tSNEGroups)) {
                for(i in seq(bc$sampleList)) {
                    groupCols[[j]][[i]] <- checkboxInput(
                        inputId=paste0("group", j, "samp", i), label="",
                        value=reactTSSamp$d[[j]][[i]]
                    )
                }
            }
        } else {
            for(j in seq(input$tSNEGroups)) {
                for(i in seq(bc$sampleList)) {
                    groupCols[[j]][[i]] <- checkboxInput(
                        inputId=paste0("group", j, "samp", i), label=""
                    )
                }
            }
        }
        checkBoxCols <- list()
        for(i in seq(input$tSNEGroups)) {
            if(!is.null(reactGroupNames$d)
               && length(reactGroupNames$d) == input$tSNEGroups) {
                textGroup <- textInput(paste0("textGroup", i), label="",
                                       value=reactGroupNames$d[[i]])
            } else {
                textGroup <- textInput(paste0("textGroup", i), label="")
            }
            butts <- list(strong(paste("Group", i)), textGroup, groupCols[[i]])
            checkBoxCols[[i]] <- column(width=1, align="center",
                                        style="margin-top: 10px;", butts)
        }
        showModal(modalDialog(
            strong("Select samples by group for concatenation."),
            br(),
            fluidRow(
                column(width=5, align="right",
                       list(br(), br(), br(), br(), sampleCol)),
                checkBoxCols
            ),
            footer=list(
                actionButton(inputId="tsneoksampbutton", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="l"))
    })

    observeEvent(input$tsneoksampbutton, {
        reactTSSamp$d <- NULL
        reactGroupNames$d <- NULL
        for(i in seq(input$tSNEGroups)) {
            reactGroupNames$d[[i]] <- input[[paste0("textGroup", i)]]
            reactTSSamp$d[[i]] <- list()
        }
        for(j in seq(input$tSNEGroups)) {
            for(i in seq(bc$sampleList)) {
                index <- paste0("group", j, "samp", i)
                reactTSSamp$d[[j]][[i]] <- input[[index]]
            }
        }
        bc$concatSamples <- list()
        for(j in seq(input$tSNEGroups)) {
            bc$concatSamples[[j]] <- list()
        }
        for(j in seq(input$tSNEGroups)) {
            bc$concatSamples[[j]] <- which(reactTSSamp$d[[j]] == TRUE)
        }
        if(length(bc$concatSamples[[1]]) == 0) {
            reactTSSamp$d <- NULL
            reactGroupNames$d <- NULL
        } else {
            for(i in seq(bc$concatSamples)) {
                bc$concatSamples[[i]] <- sapply(
                    bc$concatSamples[[i]], function(x)
                        which(bc$fileList == bc$sampleList[x]))
            }
        }
        lens <- lapply(bc$concatSamples, function(x)
            length(x) == 0)
        if(length(grep(TRUE, unlist(lens))) > 0) {
            bc$concatSamples <- list()
            alert("Select at least one sample per group.")
        } else {
            if(is.null(reactGroupNames$d)) {
                bc$concatSamples <- list()
                alert("Please type the name of each group in the
                corresponding text boxes.")
            } else {
                names <- sapply(reactGroupNames$d, function(x) x != "")
                len <- length(reactGroupNames$d)
                if(length(grep(TRUE, names)) < length(bc$concatSamples)) {
                    bc$concatSamples <- list()
                    alert("Please type the name of each group in the
                    corresponding text boxes.")
                } else if(len > 1
                          && length(unique(reactGroupNames$d)) < len) {
                    alert("Group names cannot be identical.")
                } else {
                    removeModal()
                }
            }
        }
    })

    observeEvent(input$tSNEParameters, {
        updateCheckboxGroupInput(inputId="tsneParCheckBox",
                                 selected= isolate(reactTSCh$d))
        choices <- unlist(bc$ch)[unlist(bc$ch) != "Time"]
        showModal(modalDialog(
            strong("Select at least 2 parameters to include in the analysis"),
            checkboxGroupInput("tsneParCheckBox", label="", choices=choices,
                               width="100%"),
            footer=list(
                actionButton(inputId="tsneOkParButt", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="m"))
    })

    observeEvent(input$tsneOkParButt, {
        if(length(input$tsneParCheckBox) < 2) {
            alert("Please select at least 2 parameters.")
        } else {
            reactTSCh$d <- input$tsneParCheckBox
            removeModal()
        }
    })

    observeEvent(c(reactTSSamp$d, input$tSPar, reactTSCh$d, input$tSEvs,
                   input$tabs), {
                       if(!is.null(reactTSSamp$d)
                          && input$tSPar != ""
                          && !is.null(reactTSCh$d)
                          && input$tSEvs > 0) {
                           enable("tSNEGenerate")
                       } else {
                           disable("tSNEGenerate")
                       }
                   }
    )

    observeEvent(input$tSNEGenerate, {
        if(is.null(bc$entiretSNE)) {
            events <- as.numeric(input$tSEvs)
            verifier <- sapply(unlist(bc$concatSamples), function(x)
                nrow(gs_pop_get_data(bc$gs[x], input$tSPar)[[1]]) >= events)
            if(all(verifier)) {
                len <- length(unlist(bc$concatSamples))
                round <- round(len*events/1000*2.5/60, 1)
                minutes <- as.numeric(format(round, nsmall=1))
                showModal(modalDialog(
                    "The t-SNE tool might require longer processing time
                    than other tools in this software.",
                    br(),
                    br(),
                    "Before proceeding, make sure that your CPU is not
                    being used by heavy processes.",
                    br(),
                    br(),
                    "Concatenation of the selected samples/number of
                    events and t-SNE generation are estimated to take ",
                    strong(minutes, "-", minutes*2, " minutes."),
                    br(),
                    br(),
                    textOutput("question"),
                    footer=list(
                        actionButton(inputId="tSNEGenOk",
                                     label="Generate t-SNE",
                                     style=bc$blackB),
                        actionButton(inputId="cancelModal",
                                     label="Cancel",
                                     style=bc$whiteB)),
                    easyClose=FALSE,
                    size="m"))
            } else {
                nums <- sapply(unlist(bc$concatSamples), function(x)
                    nrow(gs_pop_get_data(bc$gs[x], input$tSPar)[[1]]))
                lowerNumber <- min(nums)
                if(lowerNumber < 1000) {
                    alert("There are not enough events (<1K) in at least
                    one of the selected samples/parent.
                    Try to change the samples or the Parent.")
                } else {
                    if(lowerNumber < 5000) {
                        maxEvents <- "1K"
                    } else if(lowerNumber < 10000) {
                        maxEvents <- "5K"
                    } else if(lowerNumber < 25000) {
                        maxEvents <- "10K"
                    } else if(lowerNumber < 50000) {
                        maxEvents <- "25K"
                    }
                    alert(paste("At least one selected samples/parent
                    don't have enough events. For the selected
                                samples/parent, only", maxEvents,
                                "events are possible."))
                }
            }
        } else {
            showModal(modalDialog(
                "Do you want to make a new t-SNE?",
                footer=list(
                    actionButton(inputId="tsnemakenew", label="Make new t-SNE",
                                 style=bc$blackB),
                    actionButton(inputId="cancelModal", label="Cancel",
                                 style=bc$whiteB)
                ),
                easyClose=TRUE,
                size="m")
            )
        }
    })

    output$question <- renderText("Do you want to proceed?")

    observeEvent(input$tsnemakenew, {
        bc$entiretSNE <- NULL
        hidetSNE()
        disable("tSNEGenerate")
        tSPlotActiv$d <- tSPlotActiv$d + 1
        removeModal()
    })

    ##main----
    observeEvent(input$tSNEGenOk, {
        setwd(bc$filePath)
        reactTSPar$d <- isolate(input$tSPar)
        pops <- gh_pop_get_descendants(bc$gs[[1]], isolate(input$tSPar),
                                       path=1)
        for(i in seq(length(unlist(bc$concatSamples)))) {
            index <- unlist(bc$concatSamples)[i]
            cs <- gs_pop_get_data(bc$gs[index], "root")
            ff <- cytoframe_to_flowFrame(cs[[1]])@exprs
            if(i == 1) {
                preConcat <- as.data.frame.array(ff)
                for(k in c(1, seq(pops) + 1)) {
                    if(k < length(c(1, seq(pops) + 1))) {
                        IDs <- gh_pop_get_indices(bc$gs[[index]],
                                                  pops[k])
                        preConcat <- cbind2(preConcat, IDs)
                    } else {
                        IDs <- gh_pop_get_indices(bc$gs[[index]],
                                                  isolate(input$tSPar))
                        preConcat <- cbind2(preConcat, IDs)
                    }
                }
                pre <- which(preConcat[ncol(preConcat)][[1]] == TRUE)
                preConcat <- preConcat[pre,]
                preConcat <- preConcat[,-length(preConcat)]
                pre <- sample.int(nrow(preConcat),
                                  as.numeric(isolate(input$tSEvs)))
                preConcat <- preConcat[pre,]
            } else {
                tempConcat <- as.data.frame.array(ff)
                for(k in c(1, seq(pops) + 1)) {
                    if(k < length(c(1, seq(pops) + 1))) {
                        IDs <- gh_pop_get_indices(bc$gs[[index]],
                                                  pops[k])
                        tempConcat <- cbind2(tempConcat, IDs)
                    } else {
                        IDs <- gh_pop_get_indices(bc$gs[[index]],
                                                  isolate(input$tSPar))
                        tempConcat <- cbind2(tempConcat, IDs)
                    }
                }
                temp <- which(tempConcat[ncol(tempConcat)][[1]] == TRUE)
                tempConcat <- tempConcat[temp,]
                tempConcat <- tempConcat[,-length(tempConcat)]
                temp <- sample.int(nrow(tempConcat),
                                   as.numeric(isolate(input$tSEvs)))
                tempConcat <- tempConcat[temp,]
                preConcat <- rbind2(preConcat, tempConcat)
            }
        }
        for(i in seq(bc$fluoCh)) {
            trans <- flowjo_biexp(channelRange=4096, maxValue=262144,
                                  pos=bc$customAxis[[i]][1],
                                  neg=bc$customAxis[[i]][2],
                                  widthBasis=bc$customAxis[[i]][3],
                                  inverse=TRUE)
            preConcat[,bc$fluoCh[i]] <- trans(preConcat[,bc$fluoCh[i]])
        }
        bc$concat <- NULL
        for(i in seq(isolate(input$tsneParCheckBox))) {
            if(i == 1) {
                bc$concat <- preConcat[,isolate(input$tsneParCheckBox)[1]]
            } else {
                index <- isolate(input$tsneParCheckBox)[i]
                bc$concat <- cbind2(bc$concat, preConcat[,index])
            }
        }
        for(i in rev(ncol(preConcat) - (seq(pops) - 1))) {
            bc$concat <- cbind2(bc$concat, preConcat[,i])
        }
        bc$concat <- as.data.frame(bc$concat)
        colnames(bc$concat) <- c(isolate(input$tsneParCheckBox), pops)
        hide("tSNEGenOk")
        hide("cancelModal")
        hide("question")
        withProgress(message="Generating t-SNE...", detail="", value=0,
                     max=100, {
                         len <- length(isolate(input$tsneParCheckBox))
                         index <- seq_len(len)
                         df <- as.data.frame(bc$concat[,index])
                         bc$entiretSNE <- Rtsne(df, check_duplicates=FALSE,
                                                verbose=TRUE)
                         bc$entiretSNE <- bc$entiretSNE$Y
                     })
        showtSNE()
        removeModal()
        bc$tSNEListofGroups <- list()
        for(i in seq(bc$concatSamples)) {
            bc$tSNEListofGroups[[i]] <- i
            names(bc$tSNEListofGroups)[[i]] <- reactGroupNames$d[[i]]
        }
        bc$tSNEListofSamples <- list()
        for(i in seq(length(unlist(bc$concatSamples)))) {
            bc$tSNEListofSamples[[i]] <- i
            index <- unlist(bc$concatSamples)[i]
            names(bc$tSNEListofSamples)[[i]] <- bc$fileList[index]
        }
        if(length(bc$tSNEListofGroups) < 2
           && length(bc$tSNEListofSamples) < 2) {
            choices <- c("Heatmap", "Overlay Populations")
            updateSelectInput(inputId="tSNEMode", choices=choices)
        } else {
            choices <- c("Heatmap", "Overlay Groups or Samples",
                         "Overlay Populations")
            updateSelectInput(inputId="tSNEMode", choices=choices)
        }
        tempChannels <- unlist(bc$ch)[unlist(bc$ch) != "Time"]
        bc$availableparameters <- list()
        for(i in seq(input$tsneParCheckBox)) {
            inp <- input$tsneParCheckBox[i]
            bc$availableparameters[[i]] <- which(tempChannels == inp)
        }
        bc$availableparameters <- tempChannels[unlist(bc$availableparameters)]
        updateSelectInput(inputId="tSHighl", choices=
                              c(bc$availableparameters))
        updateSelectInput(inputId="tSNEMode", selected="Heatmap")
        updateSelectInput(inputId="tSGroupOrSamp", selected="All")
        updateSelectInput(inputId="tSGroupOrSampID", selected="")
    })

    observeEvent(input$tSGroupOrSamp, {
        if(exists("concatSamples", envir = bc)) {
            if(input$tSNEMode == "Heatmap"
               || input$tSNEMode == "Overlay Populations") {
                if(input$tSGroupOrSamp == "All") {
                    disable("tSGroupOrSampID")
                    updateSelectInput(inputId="tSGroupOrSampID",
                                      choices=c(""), selected="")
                } else {
                    enable("tSGroupOrSampID")
                    if(input$tSGroupOrSamp == "Group") {

                        updateSelectInput(inputId="tSGroupOrSampID",
                                          choices=bc$tSNEListofGroups,
                                          selected=1)
                        if(isolate(input$tSGroupOrSampID) == 1) {
                            tSPlotActiv$d <- tSPlotActiv$d + 1
                        }
                    } else if(input$tSGroupOrSamp == "Sample") {
                        updateSelectInput(inputId="tSGroupOrSampID",
                                          choices=bc$tSNEListofSamples,
                                          selected=1)
                        if(isolate(input$tSGroupOrSampID) == 1) {
                            tSPlotActiv$d <- tSPlotActiv$d + 1
                        }
                    }
                }
            } else if(input$tSNEMode == "Overlay Groups or Samples") {
                if(input$tSGroupOrSamp != "All") {
                    enable("tSGroupOrSampIDs")
                } else {
                    disable("tSGroupOrSampIDs")
                }
                reacttSIDs$d <- NULL
                tSPlotActiv$d <- tSPlotActiv$d + 1
            }
        }
    })

    observeEvent(input$tSNEMode, {
        if(!is.null(bc$entiretSNE)) {
            if(input$tSNEMode == "Heatmap"
               || input$tSNEMode == "Overlay Populations") {
                if(input$tSGroupOrSamp != "All") {
                    if(input$tSGroupOrSamp == "Group") {
                        choices <- bc$tSNEListofGroups
                    } else if(input$tSGroupOrSamp == "Sample") {
                        choices <- bc$tSNEListofSamples
                    }
                    if(input$tSGroupOrSampID != 1) {
                        updateSelectInput(inputId="tSGroupOrSampID",
                                          choices=choices, selected=1)
                    } else {
                        tSPlotActiv$d <- tSPlotActiv$d + 1
                    }
                    enable("tSGroupOrSampID")
                } else {
                    if(input$tSGroupOrSampID != "") {
                        updateSelectInput(inputId="tSGroupOrSampID",
                                          choices=c(""), selected="")
                    } else {
                        tSPlotActiv$d <- tSPlotActiv$d + 1
                    }
                    disable("tSGroupOrSampID") }
            } else if(input$tSNEMode == "Overlay Groups or Samples") {
                if(input$tSGroupOrSamp != "All") {
                    enable("tSGroupOrSampIDs")
                } else {
                    disable("tSGroupOrSampIDs")
                }
            }
            if(input$tSNEMode == "Heatmap") {
                hide("tSNEPopulations")
                hide("tSGroupOrSampIDs")
                shinyjs::show("tSGroupOrSampID")
                shinyjs::show("tSGroupOrSamp")
                shinyjs::show("tSHighl")
            } else if(input$tSNEMode == "Overlay Groups or Samples") {
                hide("tSHighl")
                hide("tSGroupOrSampID")
                hide("tSNEPopulations")
                shinyjs::show("tSGroupOrSampIDs")
                reacttSIDs$d <- NULL
                tSPlotActiv$d <- tSPlotActiv$d + 1
            } else if(input$tSNEMode == "Overlay Populations") {
                hide("tSHighl")
                hide("tSGroupOrSampIDs")
                shinyjs::show("tSGroupOrSampID")
                shinyjs::show("tSNEPopulations")
            }
        }
    })

    observeEvent(input$tSGroupOrSampIDs, {
        if(is.null(reacttSIDs$d)) {
            if(input$tSGroupOrSamp == "Group") {
                listofgrouporsamples <- bc$tSNEListofGroups
            } else if(input$tSGroupOrSamp == "Sample") {
                listofgrouporsamples <- bc$tSNEListofSamples
            }
        } else {
            if(input$tSGroupOrSamp == "Group") {
                listofgrouporsamples <- bc$tSNEListofGroups
            } else if(input$tSGroupOrSamp == "Sample") {
                listofgrouporsamples <- bc$tSNEListofSamples
            }
        }
        showModal(modalDialog(
            "Select 2 groups/samples",
            checkboxGroupInput("tsneIDs", label="",
                               choices=listofgrouporsamples,
                               selected=reacttSIDs$d),
            footer=list(
                actionButton(inputId="tSNESetIDs", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="m"))
        if(length(input$tsneIDs) == 0) {
            disable("tSNESetIDs")
        } else {
            if(length(input$tsneIDs) != 2) {
                disable("tSNESetIDs")
            } else {
                enable("tSNESetIDs")
            }
        }
    })

    observeEvent(input$tsneIDs, {
        if(input$tsneIDs[1] == "") {
            disable("tSNESetIDs")
        } else {
            if(length(input$tsneIDs) != 2) {
                disable("tSNESetIDs")
            } else {
                enable("tSNESetIDs")
            }
        }
    })

    observeEvent(input$tSNESetIDs, {
        reacttSIDs$d <- input$tsneIDs
        removeModal()
    })

    observeEvent(input$tSNEPopulations, {
        popOptions <- gh_pop_get_descendants(bc$gs[[1]], reactTSPar$d, path=1)
        showModal(modalDialog(
            "Select at least 2-8 populations.",
            checkboxGroupInput("tSNEPopIDs", label="", choices=popOptions,
                               selected=reacttSNEPops$d),
            footer=list(
                actionButton(inputId="tSNESetPopIDs", label="Ok",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="m"))
        if(length(input$tSNEPopIDs) == 0) {
            disable("tSNESetPopIDs")
        } else {
            if(length(input$tSNEPopIDs) < 2 && length(input$tSNEPopIDs) > 8) {
                disable("tSNESetPopIDs")
            } else {
                enable("tSNESetPopIDs")
            }
        }
    })

    observeEvent(input$tSNEPopIDs, {
        if(input$tSNEPopIDs[1] == "") {
            disable("tSNESetPopIDs")
        } else {
            if(length(input$tSNEPopIDs) < 2 && length(input$tSNEPopIDs) > 8) {
                disable("tSNESetPopIDs")
            } else {
                enable("tSNESetPopIDs")
            }
        }
    })

    observeEvent(input$tSNESetPopIDs, {
        reacttSNEPops$d <- input$tSNEPopIDs
        removeModal()
    })

    output$tSNEPlot <- renderPlot({
        tSPlotActiv$d
        input$tSHighl
        input$tSGroupOrSampID
        reacttSIDs$d
        if(!is.null(bc$entiretSNE)) {
            tSNEGen(input$tSNEDotSize, isolate(input$tSNEMode),
                    reacttSNEPops$d, isolate(input$tSGroupOrSamp),
                    isolate(input$tSGroupOrSampID), input$tSHighl,
                    reacttSIDs$d, isolate(input$tsneIDs), reacttSNEPops$d)


            output$showingtSNEEvents <- renderText({
                paste0("Plotted events: ", format(bc$toFormat, big.mark=","))
            })
            enable("savetSNEPlot")
            enable("tSNEGenerate")
        } else {
            disable("savetSNEPlot")
        }
    }, height=440, width=800)

    output$savetSNEPlot <- downloadHandler(
        "t-SNE plot.png",
        function(file) {
            png(file, units="in", height=6, width=11.02, res=300)
            tSNEGen(input$tSNEDotSize, isolate(input$tSNEMode),
                    reacttSNEPops$d, isolate(input$tSGroupOrSamp),
                    isolate(input$tSGroupOrSampID), input$tSHighl,
                    reacttSIDs$d, isolate(input$tsneIDs), reacttSNEPops$d)
            dev.off()
        })

    #Results----
    ##left----
    observeEvent(input$rsStat, {
        if(isolate(input$rsStat) == "Freq. of parent"
           || isolate(input$rsStat) == "Freq. of total"
           || isolate(input$rsStat) == "Count") {
            disable("rsParent")
            updateSelectInput(inputId="rsCh",
                              choices=c("", bc$completeFluoCh))
            disable("rsCh")
        } else {
            if(isolate(input$rsStat) == "Freq. of...") {
                updateSelectInput(inputId="rsCh",
                                  choices=c("", bc$completeFluoCh))
                disable("rsCh")
            } else {
                updateSelectInput(inputId="rsCh",
                                  choices=c("", bc$completeFluoCh))
                enable("rsCh")
            }
        }
    })

    observeEvent(c(input$rsStat, input$rsParent, input$rsPop, input$rsCh), {
        if(isolate(input$rsPop) != "") {
            if(isolate(input$rsStat) == "Freq. of parent"
               || isolate(input$rsStat) == "Freq. of total"
               || isolate(input$rsStat) == "Count") {
                enable("addResult")
            } else {
                if(isolate(input$rsParent) != ""
                   || isolate(input$rsCh) != "") {
                    enable("addResult")
                } else {
                    disable("addResult")
                }
            }
        } else {
            disable("addResult")
        }
    })

    observeEvent(c(input$rsStat, input$rsPop), {
        if(!is.null(bc$gs)) {
            index <- which(gs_get_pop_paths(bc$gs, path=1) == input$rsPop)
            fullPath <- gs_get_pop_paths(bc$gs)[index]
            possibleParents <- strsplit(fullPath, split="/")
            if(length(possibleParents) > 0) {
                index <- seq_len(length(possibleParents[[1]]) - 1)
                possibleParents <- possibleParents[[1]][index]
            }
            if(isolate(input$rsStat) == "Freq. of..."
               && isolate(input$rsPop) != "") {
                updateSelectInput(inputId="rsParent",
                                  choices=c(possibleParents))
                enable("rsParent")
            } else {
                updateSelectInput(inputId="rsParent", selected="")
                disable("rsParent")
            }
        }
    })

    observeEvent(input$addResult, {
        updateSelectInput(inputId="rsCh", choices=c("", bc$completeFluoCh))
        if(is.na(bc$results[[1]][1])) {
            resultID <- 1
        } else {
            resultID <- ncol(bc$results) + 1
        }
        popPaths <- gs_get_pop_paths(bc$gs)
        popShortPaths <- gs_get_pop_paths(bc$gs, path=1)
        subpop <- popPaths[which(popShortPaths == input$rsPop)]
        popStats <- gs_pop_get_count_fast(bc$gs, subpopulations=subpop)
        id <- bc$onlySampleIDs
        if(input$rsStat == "Freq. of parent") {
            num <- popStats[[4]]*100/popStats[[5]]
            bc$results[[resultID]] <- sprintf("%.1f", num)[id]
        } else if(input$rsStat == "Freq. of total") {
            nums <- gs_pop_get_count_fast(bc$gs, "freq",
                                          subpopulations=subpop)
            bc$results[[resultID]] <- sprintf("%.1f",
                                              nums$Count$Frequency*100)[id]
        } else if(input$rsStat == "Count") {
            bc$results[[resultID]] <- popStats[[4]][id]
        } else if(input$rsStat == "Freq. of...") {
            popPaths <- gs_get_pop_paths(bc$gs)
            popShortPaths <- gs_get_pop_paths(bc$gs, path=1)
            subpop <- popPaths[which(popShortPaths == input$rsParent)]
            parentstats <- gs_pop_get_count_fast(bc$gs, subpopulations=subpop)
            num <- popStats[[4]]*100/parentstats[[4]]
            bc$results[[resultID]] <- sprintf("%.1f", num)[id]
        } else if(input$rsStat == "Median") {
            index <- which(input$rsCh == bc$fluoCh) + 2
            stats <- gs_pop_get_stats(bc$gs, input$rsPop, pop.MFI)
            bc$results[[resultID]] <- sprintf("%.1f", stats[[index]])
        }
        if(input$rsStat == "Median") {
            name1 <- "MFI"
        } else {
            name1 <- input$rsStat }
        if(input$rsParent != "") {
            name2 <- paste0(input$rsParent)
        } else {
            name2 <- ""
        }
        if(input$rsCh != "") {
            name4 <- paste0(" (", input$rsCh, ")")
        } else {
            name4 <- ""
        }
        colnames(bc$results)[resultID] <- paste0(name1, name2, "/",
                                                 input$rsPop, name4)
    })

    ##main----
    output$results <- renderTable(rownames=TRUE, spacing="xs", {
        input$addResult
        input$deleteModal
        input$tabs
        if(!is.na(bc$results[,1][[1]])) {
            enable("exportTable")
            enable("editResult")
        }
        allPops <- gs_get_pop_paths(bc$gs, path=1)
        allPops <- allPops[allPops != allPops[1]]
        updateSelectInput(inputId="rsPop", choices=c("", allPops))
        bc$results
    })

    observeEvent(input$editResult, {
        choices <- setNames(as.list(seq(ncol(bc$results))),
                            colnames(bc$results))
        showModal(modalDialog(
            "Select the rows to be deleted.",
            checkboxGroupInput("checkdelete", label="", choices=choices),
            footer=list(
                actionButton(inputId="deleteModal", label="Delete",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=TRUE,
            size="s")
        )
    })

    observeEvent(input$deleteModal, {
        bc$results <- bc$results[-as.integer(input$checkdelete)]
        removeModal()
        if(ncol(bc$results) == 0) {
            bc$results[1] <- NA
            colnames(bc$results) <- NA
            disable("exportTable")
            disable("editResult")
        }
    })

    output$exportTable <- downloadHandler(
        "Results.csv",
        function(file) {
            write.csv(bc$results, file)
        })

    #Between Tabs----
    observeEvent(input$tabs, {
        if(input$tabs == "compTab") {
            bc$uncompGSfortable <- gs_clone(bc$uncompGS)
            for(i in seq(bc$fluoCh)) {
                trans <- flowjo_biexp_trans(channelRange=4096, maxValue=262144,
                                            pos=bc$customAxis[[i]][1],
                                            neg=bc$customAxis[[i]][2],
                                            widthBasis=bc$customAxis[[i]][3])
                bc$uncompGSfortable <- transform(bc$uncompGSfortable,
                                                 transformerList(bc$fluoCh[i],
                                                                 trans))
            }
            hide("saveMatrix")
            hide("cancelMatrix")
            enable("createMatrix")
            reactReadOnly$d <- TRUE
            cytomDef <- bc$compDFs["Cytometer-defined"]
            bc$hTable$d <- as.data.frame(cytomDef[[1]]*100)
            output$comp <- renderRHandsontable({
                mat <- bc$compDFs[input$previewMatrix][[1]]*100
                visualtable <- rhandsontable(mat, rowHeaderWidth=100,
                                             stretchH="all",
                                             readOnly=reactReadOnly$d)%>%
                    hot_validate_numeric(cols=seq(bc$fluoCh))%>%
                    hot_context_menu(allowRowEdit=FALSE, allowColEdit=FALSE)%>%
                    hot_cols(format="0", renderer="
            function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.NumericRenderer.apply(this,
                        arguments);
                     heatmap = ['#FFFFCC66', '#FFEDA066', '#FED97666',
                    '#FEB24C66', '#FD8D3C66', '#FC4E2A66', '#E31A1C66',
                    '#BD002666', '#80002666']
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
                      td.style.background = 'white'
                    }
            }
                   ")
                for(i in seq(bc$fluoCh)) {
                    cell <- list(row=i - 1, col=i - 1, readOnly=TRUE)
                    visualtable$x$cell <- c(visualtable$x$cell, list(cell))
                }
                onRender(visualtable, changeHook)
            })
        } else if(input$tabs == "prolifTab") {
            popPaths <- gs_get_pop_paths(bc$gs, path=1)
            if("undivided" %in% popPaths) {
                hide("prolifStart")
                shinyjs::show("prolifPlot")
            } else {
                hide("prolifPlot")
                shinyjs::show("prolifStart")
            }
        } else if(input$tabs == "resultTab") {
            updateSelectInput(inputId="rsStat", selected="")
            updateSelectInput(inputId="rsParent", selected="")
            updateSelectInput(inputId="rsPop", selected="")
            updateSelectInput(inputId="rsCh", selected="")
            disable("rsParent")
            disable("rsCh")
        }
    })

    #File----
    shinyFileChoose(input, "BCytoFileLoad", roots=c(Home= path_home()),
                    session=session, restrictions=system.file(package="base"),
                    filetypes="RData")
    shinyDirChoose(input, "directory", roots=c(Home= path_home()),
                   session=session, restrictions=system.file(package="base"),
                   allowDirCreate=FALSE, filetypes="fcs")

    observeEvent(input$directory, {
        if(!is.integer(input$directory[1])) {
            reactDir$d <- input$directory
        }
    })

    observeEvent(reactDir$d, {
        withProgress(message="Please wait...", detail="", value=0, max=100, {
            pars <- c(Home= path_home(), "R Installation"=R.home(),
                      getVolumes())
            bc$tempFilePath <- parseDirPath(pars, reactDir$d)
            bc$tempFileList <- list.files(bc$tempFilePath, pattern=".fcs")
            if(length(bc$tempFileList) == 0) {
                alert("No FCS file was found in the selected folder.")
            } else {
                setwd(bc$tempFilePath)
                validityOut <- vapply(seq(bc$tempFileList), function(x){
                        list(tryCatch({
                            keyword(load_cytoframe_from_fcs(bc$tempFileList[x])
                            )
                        }))
                    }, list(1))
                verifier <- vapply(validityOut, function(x)
                    length(x) > 0, TRUE)
                if(all(verifier)) {
                    channelLengths <- sapply(bc$tempFileList, function(x)
                        length(colnames(load_cytoframe_from_fcs(x))))
                    if(length(unique(channelLengths)) == 1) {
                        tubeNameList <- vapply(validityOut, function(x)
                            list(x$`TUBE NAME`), list(1))
                        verifier <- vapply(tubeNameList, function(x)
                            !is.null(x), TRUE)
                        if(all(verifier)) {
                            bc$tempCompIDs <- list()
                            for(i in seq(bc$tempFileList)) {
                                index <- length(bc$tempCompIDs) + 1
                                nam <- tubeNameList[[i]]
                                tubeNam <- strsplit(nam, " ")[[1]]
                                if(length(grep("Control", tubeNam)) > 0) {
                                    sta <- grep("Stained", tubeNam)
                                    unSta <- grep("Unstained", tubeNam)
                                    if(length(sta) > 0
                                       || length(unSta) > 0) {
                                        bc$tempCompIDs[[index]] <- i
                                    }
                                }
                            }
                            bc$tempCompIDs <- unlist(bc$tempCompIDs)
                            bc$tempCompIDs <- bc$tempCompIDs
                            listSeq <- seq(bc$tempFileList)
                            bc$finalSelect <- list(as.list(listSeq),
                                                   as.list(listSeq))
                            id <- bc$tempCompIDs
                            if(!is.null(id)) {
                                bc$tempOnlySampIDs <- listSeq[-id]
                                bc$finalSelect[[1]][id] <- TRUE
                                bc$finalSelect[[1]][-id] <- FALSE
                                bc$finalSelect[[2]][listSeq[-id]] <- TRUE
                                bc$finalSelect[[2]][-listSeq[-id]] <- FALSE
                                if(length(listSeq[-id]) == 0) {
                                    bc$finalSelect[[2]][id] <- FALSE
                                }
                            } else {
                                bc$finalSelect[[1]] <- lapply(
                                    bc$finalSelect[[1]], function(x) FALSE)
                                bc$finalSelect[[2]] <- lapply(
                                    bc$finalSelect[[2]], function(x) TRUE)
                            }
                            sampleColum <- list()
                            for(i in seq(bc$tempFileList)) {
                                sampleColum[[i]] <- list(br(), br(),
                                                         bc$tempFileList[i])
                            }
                            groupColums <- list(list(), list())
                            for(j in seq(2)) {
                                for(i in seq(bc$tempFileList)) {
                                    groupColums[[j]][[i]] <- checkboxInput(
                                        inputId=paste0("column", j, "samp", i),
                                        label="",
                                        value=bc$finalSelect[[j]][[i]]
                                    )
                                }
                            }
                            checkBoxColumns <- list(
                                column(
                                    width=4, align="center",
                                    style="margin-top: 10px;",
                                    list(strong("Compensation Controls"),
                                         groupColums[[1]])))
                            checkBoxColumns[[2]] <- column(
                                width=2, align="center",
                                style="margin-top: 10px;",
                                list(strong("Samples"), groupColums[[2]]))
                            showModal(modalDialog(
                                strong("Please indicate/confirm compensation
                                   controls and sample files."),
                                br(),
                                fluidRow(
                                    column(
                                        width=6, align="right",
                                        list(sampleColum)),
                                    checkBoxColumns),
                                footer=list(
                                    actionButton(inputId="okSelectCompAndSamp",
                                                 label="Ok",
                                                 style=bc$blackB),
                                    actionButton(inputId="cancelModal",
                                                 label="Cancel",
                                                 style=bc$whiteB)),
                                easyClose=FALSE,
                                size="l"))
                        } else {
                            tubeList <- vapply(validityOut, function(x)
                                list(x$`TUBE NAME`), list(1))
                            verifier <- vapply(tubeList, function(x)
                                !is.null(x), TRUE)
                            if(all(tubeList)) {
                                alert("Loading failed... At least one FCS
                                file had its integrity altered. BCyto
                                supports only original unmodified FCS files.")
                            } else {
                                alert("Loading failed... At least one FCS
                                file is of an unsupported version
                                      or is corrupted.")
                            }
                        }
                    } else {
                        alert("Loading failed... This software requires all
                        FCS files to be consistent with each other in the
                              number of channels.")
                    }
                } else {
                    alert("Loading failed... At least one FCS file seems
                          to be corrupted.")
                }
            }
            reactDir$d <- NULL
        })
    })

    observeEvent(input$okSelectCompAndSamp, {
        for(j in seq(2)) {
            for(i in seq(bc$tempFileList)) {
                index <- paste0("column", j, "samp", i)
                bc$finalSelect[[j]][[i]] <- input[[index]]
            }
        }
        displayAlert <- NULL
        for(i in seq(bc$tempFileList)) {
            if(bc$finalSelect[[1]][[i]] == FALSE
               && bc$finalSelect[[2]][[i]] == FALSE) {
                displayAlert <- "At least one file has not been checked."
            } else if(bc$finalSelect[[1]][[i]] == TRUE
                      && bc$finalSelect[[2]][[i]] == TRUE) {
                displayAlert <- "At least one file has been checked twice."
            }
        }
        verifier <- sapply(bc$finalSelect[[2]], function(x)
            x == FALSE)
        if(all(verifier)) {
            displayAlert <- "Only compensation control files were loaded.
            This software requires at least a sample file to be loaded."
        }
        if(is.null(displayAlert)) {
            withProgress(message="Please wait...", detail="", value=0,
                         max=100, {
                             hide("okSelectCompAndSamp")
                             hide("cancelModal")
                             bc$filePath <- bc$tempFilePath
                             bc$fileList <- bc$tempFileList
                             preList <- unlist(bc$finalSelect[[1]])
                             bc$compControlIDs <- seq(bc$fileList)[preList]
                             preList <- unlist(bc$finalSelect[[2]])
                             bc$onlySampleIDs <- seq(bc$fileList)[preList]
                             initializing()
                             initInpUpd()
                             plotActivator$d <- plotActivator$d + 1
                             hierarchyActivator$d <- hierarchyActivator$d + 1
                             matChoices()
                             updateSelectInput(inputId="previewMatrix",
                                               choices=bc$choices,
                                               selected=bc$appliedMatrix)
                             updateSelectInput(inputId="bgType", selected="")
                             updateSelectInput(inputId="ovTyp", selected="")
                             updateSelectInput(inputId="tSEvs", selected=0)
                             disableTabs()
                             hidetSNE()
                             disable("editResult")
                             disable("exportTable")
                             if(reactTestData$d == FALSE) {
                                 enable("saveFile")
                             } else {
                                 reactTestData$d <- FALSE
                             }
                             enable("tSNEGenerate")
                             bc$betweenPeaks <- NULL
                             bc$refPeaks <- NULL
                             bc$entiretSNE <- NULL
                             bc$tSNEListofGroups <- NULL
                             bc$tSNEListofSamples <- NULL
                             bc$concat <- NULL
                             bc$concatSamples <- NULL
                             bc$samples <- NULL
                             reactTSPar$d <- NULL
                             reactTSCh$d <- NULL
                             reactGroupNames$d <- NULL
                             reacttSIDs$d <- NULL
                             reacttSNEPops$d <- NULL
                             reactOvSamples$d <- NULL
                             tSPlotActiv$d <- tSPlotActiv$d + 1
                             if(length(bc$fluoCh) >= 10) {
                                 reactHeight$d <- 850
                             } else if(length(bc$fluoCh) < 10
                                       && length(bc$fluoCh) > 4) {
                                 reactHeight$d <- length(bc$fluoCh)*85
                             } else if(length(bc$fluoCh) <= 4) {
                                 reactHeight$d <- 350
                             }
                             alert("Files have been loaded successfully.")
                             removeModal()
                         })
        } else {
            if(!is.null(bc$filePath)){
                setwd(bc$filePath)
            }
            alert(displayAlert)
        }
    })

    observeEvent(input$loadTestData, {
        reactTestData$d <- TRUE
        file1 <- "Cell maturation assay - 24 mb
        (FlowRepository ID FR-FCM-Z4LZ)"
        file2 <- "Proliferation assay - 18.5 mb
        (FlowRepository ID FR-FCM-Z4LY)"
        showModal(modalDialog(
            strong("Select the data set to be downloaded."),
            br(),
            radioButtons("radioTestData", label="",
                         choiceValues=c("test-data-1", "test-data-2"),
                         choiceNames=c(file1, file2), width="100%"),
            footer=list(
                actionButton(inputId="downloadTestData", label="Download",
                             style=bc$blackB),
                actionButton(inputId="cancelModal", label="Cancel",
                             style=bc$whiteB)
            ),
            easyClose=FALSE,
            size="m"))
    })

    observeEvent(input$downloadTestData, {
        hide("downloadTestData")
        hide("cancelModal")
        if(input$radioTestData == "test-data-1") {
            bc$fileList <- c("Compensation Controls_APC Stained Control_013
                             .fcs",
                             "Compensation Controls_FITC Stained Control_011
                             .fcs",
                             "Compensation Controls_PE Stained Control_014
                             .fcs",
                             "Compensation Controls_PE-Cy7 Stained Control_015
                             .fcs",
                             "Compensation Controls_Unstained Control_010
                             .fcs",
                             "Compensation Controls_V 450 Stained Control_012
                             .fcs",
                             "Specimen_001_FMO-CD11c_003.fcs",
                             "Specimen_001_FMO-CD86_005.fcs",
                             "Specimen_001_FMO-Gr1_002.fcs",
                             "Specimen_001_FMO-MHCII_004.fcs",
                             "Specimen_001_FMO-Viability_001.fcs",
                             "Specimen_001_iDC_001_006.fcs",
                             "Specimen_001_iDC_002_007.fcs",
                             "Specimen_001_mDC_001_008.fcs",
                             "Specimen_001_mDC_002_009.fcs")
            bc$compControlIDs <- seq_len(6)
            bc$onlySampleIDs <- c(7:15)
            datasetFile <- "dataset-1.RData"
        } else {
            bc$fileList <- c("Compensation Controls_APC-Cy7 Stained Control_080
                             .fcs",
                             "Compensation Controls_FITC Stained Control_078
                             .fcs",
                             "Compensation Controls_PE Stained Control_081
                             .fcs",
                             "Compensation Controls_Unstained Control_077
                             .fcs",
                             "Compensation Controls_V 450 Stained Control_079
                             .fcs",
                             "Specimen_001_0-00ug-mL_002_030.fcs",
                             "Specimen_001_0-01ug-mL_003_034.fcs",
                             "Specimen_001_0-05ug-mL_002_036.fcs",
                             "Specimen_001_1-00ug-mL_001_044.fcs",
                             "Specimen_001_5-00ug-mL_002_048.fcs",
                             "Specimen_001_FMO-CD4_003.fcs",
                             "Specimen_001_FMO-CD44_004.fcs",
                             "Specimen_001_FMO-CFSE_001.fcs",
                             "Specimen_001_FMO-Viability_002.fcs")
            bc$compControlIDs <- seq_len(5)
            bc$onlySampleIDs <- c(6:14)
            datasetFile <- "dataset-2.RData"
        }
        bc$filePath <- tempdir()
        verifier <- vapply(bc$fileList, function(x)
            length(grep(x, list.files(bc$filePath))) > 0, TRUE)
        if(all(verifier) == FALSE) {
            adjustedFileList <- list()
            for(i in seq(bc$fileList)) {
                if(grepl(" ", bc$fileList[i])) {
                    adjustedFileList[[i]] <- gsub(" ", "%20", bc$fileList[i])
                } else {
                    adjustedFileList[[i]] <- bc$fileList[i]
                }
            }
            adjustedFileList <- unlist(adjustedFileList)
            dat <- data.frame(x=numeric(0), y=numeric(0))
            message <- "Downloading"
            max <- length(bc$fileList) + 1
            withProgress(message=message, detail="file 0", value=0, max=max, {
                for(i in seq(bc$fileList)) {
                    df <- data.frame(x=stats::rnorm(1), y=stats::rnorm(1))
                    dat <- rbind(dat, df)
                    incProgress(1, detail=paste("file", i))
                    Sys.sleep(0.1)
                    download.file(paste0("https://github.com/BonilhaCaio/",
                                         input$radioTestData, "/raw/main/",
                                         adjustedFileList[i]),
                                  destfile=file.path(bc$filePath,
                                                     bc$fileList[i]))
                }
                download.file(paste0("https://github.com/BonilhaCaio/",
                                     input$radioTestData, "/raw/main/",
                                     datasetFile),
                              destfile=file.path(bc$filePath, datasetFile))
            })
        }
        reactBCytoFile$d <- datasetFile
        if(input$radioTestData == "test-data-1") {
            js$disableTab("prolifTab")
        } else {
            js$enableTab("prolifTab") }
        js$enableTab("ancestryTab")
        js$enableTab("overlayTab")
        js$enableTab("tsneTab")
        js$enableTab("resultTab")
    })

    observeEvent(input$BCytoFileLoad, {
        if(!is.integer(input$BCytoFileLoad[1])) {
            reactBCytoFile$d <- input$BCytoFileLoad
        }
    })

    observeEvent(reactBCytoFile$d, {
        if(reactBCytoFile$d[[1]] == "dataset-1.RData"
           || reactBCytoFile$d[[1]] == "dataset-2.RData") {
            verifier <- TRUE
        } else {
            tempFilePath <- parseFilePaths(c(Home= path_home(),
                                             "R Installation"=R.home(),
                                             getVolumes()),
                                           reactBCytoFile$d)
            tempStrSplit <- strsplit(tempFilePath[[4]][[1]], "")[[1]]
            end <- length(tempStrSplit)
            start <- end - length(strsplit(tempFilePath[[1]], "")[[1]])
            bc$filePath <- paste(tempStrSplit[-start:-end], collapse="")
            load(paste0(bc$filePath, "/", tempFilePath[[1]]), tempEnv)
            verifier <- vapply(tempEnv$fileList, function(x)
                length(grep(x, list.files(bc$filePath))) > 0, TRUE)
        }
        if(all(verifier)) {
            bc$loadingFile <- TRUE
            withProgress(message="Please wait...", detail="", value=0,
                         max=100, {
                             setwd(bc$filePath)
                             setOne <- "dataset-1.RData"
                             setTwo <- "dataset-2.RData"
                             if(reactBCytoFile$d[[1]] == setOne
                                || reactBCytoFile$d[[1]] == setTwo) {
                                 load(paste0(bc$filePath, "/",
                                             reactBCytoFile$d), bc)
                             } else {
                                 load(paste0(bc$filePath, "/",
                                             tempFilePath[[1]]), bc)
                             }
                             bc$sampleList <- bc$fileList[bc$onlySampleIDs]
                             cs <- load_cytoset_from_fcs(bc$fileList)
                             bc$uncompGS <- GatingSet(cs)
                             bc$gs <- gs_clone(bc$uncompGS)
                             compensate(bc$gs,
                                        bc$compDFs[bc$appliedMatrix][[1]])
                             indices <- grepl("SC", colnames(bc$gs)) != TRUE
                             preFluoCh <- colnames(bc$gs)[indices]
                             bc$fluoCh <- preFluoCh[preFluoCh != "Time"]
                             completeFluoCh <- unlist(bc$ch)
                             ind <- completeFluoCh != "Time"
                             completeFluoCh <- completeFluoCh[ind]
                             indices <- grepl("SC", completeFluoCh) != TRUE
                             bc$completeFluoCh <- completeFluoCh[indices]
                             preD <- as.data.frame(bc$compDFs[[1]]*100)
                             bc$hTable <- reactiveValues(d=preD)
                             for(i in seq(bc$fluoCh)) {
                                 pre <- bc$customAxis[[i]]
                                 trans <- flowjo_biexp_trans(channelRange=4096,
                                                             maxValue=262144,
                                                             pos=pre[1],
                                                             neg=pre[2],
                                                             widthBasis=pre[3])
                                 tList <- transformerList(bc$fluoCh[i], trans)
                                 bc$gs <- transform(bc$gs, tList)
                             }
                             for(i in seq(bc$popGates)) {
                                 gs_pop_add(bc$gs, bc$popGates[[i]],
                                            parent=bc$popParents[i])
                             }
                             recompute(bc$gs)
                             reactShowAxis$d <- bc$plotShowAxis
                             reactAxisFont$d <- bc$plotAxisFont
                             reactShowGateName$d <- bc$plotShowGateName
                             reactGateFont$d <- bc$plotGateFont
                             matChoices()
                             updateSelectInput(inputId="previewMatrix",
                                               choices=bc$choices,
                                               selected=bc$appliedMatrix)
                             pops <- gs_get_pop_paths(bc$gs, path=1)
                             choices <- c("", bc$sampleList)
                             updateSelectInput(inputId="bgPreviewSample",
                                               choices=choices,
                                               selected=bc$bgPreviewSample)
                             updateSelectInput(inputId="bgType",
                                               selected=bc$bgType)
                             updateSelectInput(inputId="bgPop",
                                               choices=c("", pops[-1]),
                                               selected=bc$bgPop)
                             showOvPlot$d <- bc$ovBoolean
                             updateSelectInput(inputId="ovTyp",
                                               selected=bc$ovTyp)
                             updateSelectInput(inputId="ovTon",
                                               selected=bc$ovTon)
                             updateSelectInput(inputId="ovY", choices=c(bc$ch),
                                               selected=bc$ovY)
                             updateSelectInput(inputId="ovX", choices=c(bc$ch),
                                               selected=bc$ovX)
                             updateSelectInput(inputId="ovP", choices=pops,
                                               selected=bc$ovP)
                             reactOvSamples$d <- bc$ovSamples
                             reactOvShowAxis$d <- bc$ovShowAxis
                             reactOvAxisFont$d <- bc$ovAxisFont
                             prolifReady$d <- bc$prolifBool
                             updateSelectInput(inputId="prolifSButt",
                                               choices=bc$sampleList,
                                               selected=bc$prolifSButt)
                             if(bc$prolifBool == TRUE) {
                                 updateCheckboxInput(inputId="prolifLabel",
                                                     value=bc$prolifLabel)
                                 updateCheckboxInput(inputId="prolifGrid",
                                                     value=bc$prolifGrid)
                             }
                             reactTSPar$d <- bc$tSPar
                             reactTSCh$d <- bc$tSNEParameters
                             if(is.null(bc$tSEvs)) {
                                 bc$tSEvs <- 0
                             }
                             updateSelectInput(inputId="tSEvs",
                                               selected=bc$tSEvs)
                             reactGroupNames$d <- bc$tSNEGroupNames
                             if(length(bc$tSNEListofGroups) < 2
                                && length(bc$tSNEListofSamples) < 2) {
                                 choices <- c("Heatmap", "Overlay Populations")
                                 updateSelectInput(inputId="tSNEMode",
                                                   choices=choices,
                                                   selected=bc$tSNEMode)
                             } else {
                                 choices <- c("Heatmap",
                                              "Overlay Groups or Samples",
                                              "Overlay Populations")
                                 updateSelectInput(inputId="tSNEMode",
                                                   choices=choices,
                                                   selected=bc$tSNEMode)
                             }
                             updateSelectInput(inputId="tSHighl",
                                               selected=bc$tSHighl)
                             updateSelectInput(inputId="tSGroupOrSamp",
                                               selected=bc$tSGroupOrSamp)
                             updateSelectInput(inputId="tSGroupOrSampID",
                                               selected=bc$tSGroupOrSampID)
                             reacttSIDs$d <- bc$tSGroupOrSampIDs
                             reacttSNEPops$d <- bc$tSNEPops
                             updateSliderInput(inputId="tSNEDotSize",
                                               value=bc$tSNEDotSize)
                             if(!is.null(bc$entiretSNE)) {
                                 showtSNE()
                                 ind <- unlist(bc$ch) != "Time"
                                 tempChannels <- unlist(bc$ch)[ind]
                                 bc$availableparameters <- list()
                                 for(i in seq(bc$tSNEParameters)) {
                                     preId <- bc$tSNEParameters[i]
                                     index <- which(tempChannels == preId)
                                     bc$availableparameters[[i]] <- index
                                 }
                                 indices <- unlist(bc$availableparameters)
                                 preParam <- tempChannels[indices]
                                 bc$availableparameters <- preParam
                                 choices <- c(bc$availableparameters)
                                 updateSelectInput(inputId="tSHighl",
                                                   choices=choices)
                                 tSPlotActiv$d <- tSPlotActiv$d + 1
                             } else {
                                 hidetSNE()
                                 enable("tSNEGenerate")
                             }
                         })
            plotActivator$d <- plotActivator$d + 1
            bc$X <- colnames(bc$gs)[grep("FS", colnames(bc$gs))[1]]
            bc$Y <- colnames(bc$gs)[grep("SS", colnames(bc$gs))[1]]
            updateSelectInput(inputId="samp", choices=c("", bc$sampleList),
                              selected=bc$sampleList[1])
            updateSelectInput(inputId="typ", selected="Pseudocolor")
            updateSelectInput(inputId="previewSample",
                              choices=c("", bc$fileList),
                              selected=bc$sampleList[1])
            updateSelectInput(inputId="Y", choices=c(bc$ch), selected=bc$Y)
            updateSelectInput(inputId="X", choices=c(bc$ch), selected=bc$X)
            updateSelectInput(inputId="rsCh", choices=c("", bc$completeFluoCh))
            js$enableTab("plotTab")
            js$enableTab("compTab")
            disableTabs()
            removeModal()
            alert("File has been loaded successfully.")
            if(reactTestData$d == FALSE) {
                enable("saveFile")
            } else {
                reactTestData$d <- FALSE
            }
            reactPar$d <- "root"
            hierarchyActivator$d <- hierarchyActivator$d + 1
            reactBCytoFile$d <- NULL
            if(length(bc$fluoCh) >= 10) {
                reactHeight$d <- 850
            } else if(length(bc$fluoCh) < 10 && length(bc$fluoCh) > 4) {
                reactHeight$d <- length(bc$fluoCh)*85
            } else if(length(bc$fluoCh) <= 4) {
                reactHeight$d <- 350
            }
        } else {
            alert("Loading failed... Please certify this save file and all
                  FCS files connected to it are in the folder.")
        }
    })

    observeEvent(input$saveFile, {
        if(file.exists(paste0(getwd(), "/BCytoSave.RData"))) {
            showModal(modalDialog(
                paste("This will overwrite the existing BCyto file in",
                      getwd(), ". Do you want to proceed?"),
                footer=list(
                    actionButton(inputId="overwrite",
                                 label="Save",
                                 style=bc$blackB),
                    actionButton(inputId="cancelModal", label="Cancel",
                                 style=bc$whiteB)
                ),
                easyClose=FALSE,
                size="m")
            )
        } else {
            if(is.null(reactOverw$d)) {
                reactOverw$d <- 0
            } else {
                reactOverw$d <- reactOverw$d + 1
            }
        }
    })

    observeEvent(input$overwrite, {
        if(is.null(reactOverw$d)) {
            reactOverw$d <- 0
        } else {
            reactOverw$d <- reactOverw$d + 1
        }
        removeModal()
    })

    observeEvent(reactOverw$d, {
        setwd(bc$filePath)
        popPaths <- gs_get_pop_paths(bc$gs)
        if(length(popPaths) == 1) {
            popParents <- NULL
            popGates <- NULL
        } else {
            popParents <- gs_pop_get_count_fast(bc$gs[1], "freq")[,3][[1]]
            popGates <- vapply(popPaths[-1], function(x)
                list(gs_pop_get_gate(bc$gs, x)[1][[1]]),
                list(length(popPaths[-1])))
        }
        plotShowAxis <- reactShowAxis$d
        plotAxisFont <- reactAxisFont$d
        plotShowGateName <- reactShowGateName$d
        plotGateFont <- reactGateFont$d
        bgPreviewSample <- input$bgPreviewSample
        bgType <- input$bgType
        bgPop <- input$bgPop
        ovBoolean <- showOvPlot$d
        ovTyp <- input$ovTyp
        ovTon <- input$ovTon
        ovY <- input$ovY
        ovX <- input$ovX
        ovP <- input$ovP
        ovSamples <- reactOvSamples$d
        ovShowAxis <- reactOvShowAxis$d
        ovAxisFont <- reactOvAxisFont$d
        prolifBool <- prolifReady$d
        prolifSButt <- input$prolifSButt
        prolifLabel <- input$prolifLabel
        prolifGrid <- input$prolifGrid
        tSPar <- reactTSPar$d
        tSNEParameters <- reactTSCh$d
        tSEvs <- input$tSEvs
        tSNEGroupNames <- reactGroupNames$d
        tSNEMode <- input$tSNEMode
        tSHighl <- input$tSHighl
        tSGroupOrSamp <- input$tSGroupOrSamp
        tSGroupOrSampID <- input$tSGroupOrSampID
        tSGroupOrSampIDs <- reacttSIDs$d
        tSNEPops <- reacttSNEPops$d
        tSNEDotSize <- input$tSNEDotSize
        fileList <- bc$fileList
        compControlIDs <- bc$compControlIDs
        onlySampleIDs <- bc$onlySampleIDs
        ch <- bc$ch
        compDFs <- bc$compDFs
        appliedMatrix <- bc$appliedMatrix
        prolifChannel <- bc$prolifChannel
        prolifParent <- bc$prolifParent
        betweenPeaks <- bc$betweenPeaks
        refPeaks <- bc$refPeaks
        adjustedL <- bc$adjustedL
        concat <- bc$concat
        entiretSNE <- bc$entiretSNE
        concatSamples <- bc$concatSamples
        tSNEListofGroups <- bc$tSNEListofGroups
        tSNEListofSamples <- bc$tSNEListofSamples
        results <- bc$results
        customAxis <- bc$customAxis
        save(fileList, compControlIDs, onlySampleIDs, customAxis, popParents,
             popGates, plotShowAxis, plotAxisFont, plotShowGateName,
             plotGateFont, ch, #Plot
             compDFs, appliedMatrix, #Compensation
             bgPreviewSample, bgType, bgPop, #Ancestry
             ovBoolean, ovTyp, ovTon, ovY, ovX, ovP,
             ovSamples, ovShowAxis, ovAxisFont, #Overlays
             prolifChannel, prolifParent, prolifBool, betweenPeaks, refPeaks,
             adjustedL, prolifSButt, prolifLabel,
             prolifGrid, #Proliferation
             concat, entiretSNE, concatSamples, tSPar, tSNEParameters,
             tSEvs, tSNEGroupNames, tSNEMode, tSHighl,
             tSGroupOrSamp, tSGroupOrSampID, tSGroupOrSampIDs,
             tSNEPops, tSNEDotSize, tSNEListofGroups, tSNEListofSamples, #t-SNE
             results, #Results
             file=paste0(getwd(), "/BCytoSave.RData"))
        alert(paste("A BCyto file has been saved in the same folder of the
                FCS files:", bc$filePath))
    })

    observeEvent(input$fileHelp, {
        showModal(modalDialog(
            title="Help",
            "To start, load your FCS files, a BCyto file or download
            test data.",
            br(),
            br(),
            "Donwloading test data will load FCS files from
            https://github.com/BonilhaCaio/ in a temporary folder that will
            be automatically deleted when the R session is terminated.",
            br(),
            br(),
            "Saving the data will create a BCytoSave.RData file inside the
            FCS files folder.",
            br(),
            br(),
            strong("Important:"),
            "BCyto files are required to be inside the loaded FCS files
            folder before being uploaded.",
            footer=NULL,
            size="m",
            easyClose=TRUE))
    })
}
