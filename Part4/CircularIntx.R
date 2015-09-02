plotCircIntx <- function(pmatrix, threshold=1e-20, map, file) 
{
  wh <- which(pmatrix < threshold, arr.ind=T)
  wh <- cbind(wh, apply(wh, 1, function(x) pmatrix[x[1], x[2]]))
  wh <- as.data.frame(wh)
  names(wh) <- c("index1", "index2", "logP")
  wh <- wh[wh[,1]<wh[,2],]

  library(RCircos)
  png(file, units="in", res=600, width=8, height=8)

  nchr <- length(map)
  chr <- rep(paste("chr", 1:nchr, sep=""), unlist(lapply(map, function(x) length(x)-1)))
  map2 <- map
  for (i in 1:nchr) {
    map2[[i]] <- map2[[i]]+1:length(map2[[i]])/1000000
    map2[[i]] <- map2[[i]]-min(map2[[i]])
  }
  start <- unlist(lapply(map2, function(x) x[1:(length(x)-1)]))
  end <- unlist(lapply(map2, function(x) x[2:length(x)])) 
  start <- start*1000000
  end <- end*1000000

  band1 <- rep("p36.31", length(chr))
  stain1 <- rep(names(map), each=2)
  stain1 <- rep(names(map2), unlist(lapply(map2, function(x) length(x)-1)))

  circosdata <- data.frame(Chromosome=chr, ChromStart=start, ChromEnd=end, Band=band1, Stain=stain1)

  ## Base plot
  cyto.info <- circosdata
  chr.exclude <- NULL
  num.inside <- 3
  num.outside <- 0
  RCircos.Set.Core.Components(cyto.info, chr.exclude, num.inside, num.outside)
  RCircos.Set.Plot.Area()

  params <- RCircos.Get.Plot.Ideogram()
  params$BandColor <- params$ChrColor
  params$ChrColor <- rep("white", length(params$BandColor))
  params$Chromosome <- rep(paste("chr", names(map2), sep=""), unlist(lapply(map, function(x) length(x)-1)))
  RCircos.Reset.Plot.Ideogram(params)

  RCircos.Chromosome.Ideogram.Plot()
  params <- RCircos.Get.Plot.Ideogram()
  params$Chromosome <- cyto.info$Chromosome
  RCircos.Reset.Plot.Ideogram(params)

  ## Create link data frame
  pos <- unlist(map2)*1000000

  chr <- paste("chr", 1:nchr, sep="")
  chr <- rep(chr, unlist(lapply(map2, length)))
  chr1 <- chr[wh[,1]]
  chr2 <- chr[wh[,2]]
  index1 <- which(chr1!=(chr[wh[,1]+1]))
  start1 <- pos[wh[,1]]
  end1 <- pos[wh[,1]+1]
  if (length(index1)>0) {
     start1[index1] <- pos[wh[index1,1]-1]
     end1[index1] <- pos[wh[index1,1]]
  } 
  index2 <- which(chr2!=(chr[wh[,2]+1]))
  start2 <- pos[wh[,2]]
  end2 <- pos[wh[,2]+1]
  if (length(index2)>0) {
    start2[index2] <- pos[wh[index2,2]-1]
    end2[index2] <- pos[wh[index2,2]]
  }
 
  circoslink <- data.frame(Chromosome=chr1, chromStart=start1, chromEnd=end1, Chromosome.1=chr2, chromStart.1=start2, chromEnd.1=end2)

  ## set up colors - seven colors for each homology group and then red for not 
  ## within HG
  library(RColorBrewer)
#  hg <- rep(1:7, each=3)
  hg <- 1:length(map)
#  hg <- rep(hg, unlist(lapply(map2, length)))
  PlotColor <- rep("black", nrow(circoslink))
  for (i in 1:length(map))
	PlotColor[hg[wh[,1]]==hg[wh[,2]] & hg[wh[,1]]==i] <- brewer.pal(length(map), "Paired")[i]

  circoslink <-cbind(circoslink, PlotColor)

  ### To get around issues with plotting
  circoslink[circoslink[,1]=="chr1" & circoslink[,2]==0,2] <- 1500
  circoslink[circoslink[,4]=="chr1" & circoslink[,5]==0,5] <- 1500

  track.num <- 2
  RCircos.Link.Plot(circoslink, track.num, FALSE)
  ## may need to convert this to ribbon plot 

dev.off()
}

