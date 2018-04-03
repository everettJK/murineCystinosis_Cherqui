# library(rtracklayer)
# library(BSgenome.Mmusculus.UCSC.mm9)
# export(BSgenome.Mmusculus.UCSC.mm9, 'mm9.2bit')
# cmd <- 'blat mm9.2bit PMC3129560.ff PMC3129560.intSiteCaller.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead'

options(stringsAsFactors = FALSE)
b <- read.table('PMC3129560.intSiteCaller.psl', header=FALSE)
names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 
              'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')

b$percentID <- (b$matches/b$qSize)*100
b$percentCoverage <- ((b$qEnd - b$qStart)/b$qSize)*100

bf <- subset(b, percentID >= 98 & tBaseInsert <= 3)

t <- do.call(rbind, lapply(split(bf, bf$qName), function(x){
  if(nrow(x)==1){
   return(x)
  } else {
    return(data.frame())
  }
}))

library(GenomicRanges)
intSites <- makeGRangesFromDataFrame(t, seqnames.field= 'tName', start.field = 'tStart', end.field = 'tEnd', strand.field = 'strand', starts.in.df.are.0based = TRUE)

intSites <- intSites[width(intSites) >= 8]

library(dplyr)
library(hiReadsProcessor)
gr <- intSites
gr$Position <- ifelse(strand(gr) == "+", start(gr), end(gr))
gr$Break <- ifelse(strand(gr) == "+", end(gr), start(gr))
gr$Score <- 95
gr$qEnd <- width(gr)
gr.std <- clusterSites(psl.rd = gr, weight = rep(1, length(gr)),  windowSize = 5L)

start(gr.std) <- ifelse(strand(gr.std) == "+", gr.std$clusteredPosition, gr.std$Break)
end(gr.std)   <- ifelse(strand(gr.std) == "-", gr.std$clusteredPosition, gr.std$Break)
for (i in c("Position", "Break", "Score", "qEnd", "clusteredPosition", "clonecount", "clusterTopHit")) { mcols(gr.std[, i]) <- NULL }
intSites <- gr.std

intSites <- flank(intSites, -1, start = TRUE)

library(inSite)
library(parallel)
cluster <- makeCluster(30)

intSites$s <- ceiling(seq_along(intSites)/(length(intSites)/30))
intSites   <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$s), function(p){
  library(inSite)
  nearestGenomicFeature2(p, genome='mm9')
})))
intSites$s <- NULL

PMC3129560.intSiteData <- intSites
save(PMC3129560.intSiteData, file='PMC3129560.RData')







