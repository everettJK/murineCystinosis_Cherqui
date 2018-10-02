
nearestGenomicFeature <- function(query, genome='hg38', side='either', geneList=NULL){
  
  if(tolower(genome) == 'hg38'){
    subject       <- gt23::hg38.refSeqGenesGRanges
    subject.exons <- gt23::hg38.refSeqGenesGRanges.exons
  } else if (tolower(genome) == 'mm9') {
    subject       <- gt23::mm9.refSeqGenesGRanges
    subject.exons <- gt23::mm9.refSeqGenesGRanges.exons # Missing
  } else if (tolower(genome) == 'susscr3') {
    subject       <- gt23::susScr3.refSeqGenesGRanges
    subject.exons <- gt23::susScr3.refSeqGenesGRanges.exons
  } else {
    stop('There is not refSeq table for the requested genome.')
  }
  
  if(! is.null(geneList)){
    subject <- GenomicRanges::subset(subject, toupper(name2) %in% toupper(geneList))
  }
  
  # If side is not set to either, collapse the subject ranges to single positions
  if (side %in% c("5p", "3p", "midpoint")) {
    options(warn=-1)
    if (side == "5p") subject <- GenomicRanges::flank(subject, width = -1)
    if (side == "3p") subject <- GenomicRanges::flank(subject, width = -1, start = FALSE)
    if (side == "midpoint") ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
    ###subject <- subject[-GenomicRanges:::get_out_of_bound_index(subject)]
    options(warn=0)
  }
  
  options(stringsAsFactors = FALSE)
  
  query.df  <- GenomicRanges::as.data.frame(query)
  subject.df <- GenomicRanges::as.data.frame(subject)
  
  query.df$strand <- as.character(query.df$strand)
  subject.df$strand <- as.character(subject.df$strand)
  
  
  subject.exons.df <- GenomicRanges::as.data.frame(subject.exons)
  query.df$inFeature            <- FALSE
  query.df$nearestFeature       <- 'None.found'
  query.df$nearestFeatureStrand <- 'None.found'
  query.df$inFeatureExon        <- FALSE
  query.df$inFeatureSameOrt     <- FALSE
  query.df$nearestFeatureStart  <- Inf
  query.df$nearestFeatureEnd    <- Inf
  query.df$nearestFeatureDist   <- Inf
  
  o  <- suppressWarnings(GenomicRanges::nearest(query, subject, select='all', ignore.strand=TRUE))
  
  if(length(o) > 0){
    
    createCol <- function(a, b, n){
      paste0(unique(cbind(a, b))[,n], collapse=',')
    }
    
    ### browser()
    
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(
                    gene   = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 1),
                    strand = createCol(subject.df[subjectHits,]$name2, subject.df[subjectHits,][['strand']], 2),
                    hitStart = min(subject.df[subjectHits,][['start']]),
                    hitEnd   = max(subject.df[subjectHits,][['end']])) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene, strand, hitStart, hitEnd) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    query.df[a$queryHits,]$nearestFeature       <- a$gene
    query.df[a$queryHits,]$nearestFeatureStrand <- a$strand
    query.df[a$queryHits,]$nearestFeatureStart  <- a$hitStart
    query.df[a$queryHits,]$nearestFeatureEnd    <- a$hitEnd
  }
  
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    #query.df[a$queryHits,]$inFeature <- a$gene
    query.df[a$queryHits,]$inFeature <- TRUE
  }
  
  o <- suppressWarnings(GenomicRanges::distanceToNearest(query,  subject, select='all', ignore.strand=TRUE))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::top_n(-1, distance) %>%
      dplyr::ungroup() %>%
      dplyr::select(queryHits, distance) %>% 
      dplyr::distinct() %>% 
      data.frame()
    query.df[a$queryHits,]$nearestFeatureDist <- a$distance
  }
  
  query.df$nearestFeatureBoundary <- ifelse(abs(query.df$start - query.df$nearestFeatureStart) > 
                                              abs(query.df$start - query.df$nearestFeatureEnd),   
                                            query.df$nearestFeatureEnd,  
                                            query.df$nearestFeatureStart)
  
  query.df$nearestFeatureDist <- query.df$nearestFeatureDist * sign(query.df$start - query.df$nearestFeatureBoundary)
  query.df$nearestFeatureDist <- ifelse(query.df$nearestFeatureStrand=='+', query.df$nearestFeatureDist, query.df$nearestFeatureDist * -1)
  
  query.df$nearestFeatureStart    <- NULL
  query.df$nearestFeatureEnd      <- NULL
  query.df$nearestFeatureBoundary <- NULL
  
  # In exon
  o <- suppressWarnings(GenomicRanges::findOverlaps(query, subject.exons, select='all', ignore.strand=TRUE, type='any'))
  if(length(o) > 0){
    a <- dplyr::group_by(data.frame(o), queryHits) %>% 
      dplyr::mutate(gene   = paste(unique(subject.exons.df[subjectHits,]$name2), collapse=',')) %>% 
      dplyr::ungroup() %>%
      dplyr::select(queryHits, gene) %>% 
      dplyr::distinct() %>% 
      data.frame()
    
    #query.df[a$queryHits,]$inFeatureExon  <- a$gene
    query.df[a$queryHits,]$inFeatureExon  <- TRUE
  }
  
  # In TU ort
  # There may be cases where a site overlaps two or more features which have the sampe orientation.
  # ie. +,+,+ and we want to reduce these down to a single unique sign for comparison.
  a <- query.df[is.na(query.df$inFeature),]  
  b <- query.df[! is.na(query.df$inFeature),]
  
  if(nrow(a) > 0 && nrow(b) > 0){
    b$nearestFeatureStrandCmp <- unlist(lapply(strsplit(b$nearestFeatureStrand, ','), function(x){ paste(unique(x), collapse=',')}))
    b$inFeatureSameOrt <- b$strand == b$nearestFeatureStrandCmp
    b$nearestFeatureStrandCmp <- NULL
    query.df <- dplyr::bind_rows(a, b)
  }
  
  GenomicRanges::makeGRangesFromDataFrame(query.df, keep.extra.columns = TRUE)
}






createUCSCintSiteAbundTrack <- function(posid, abund, subject, title='intSites', outputFile='track.ucsc', position=NA){
  o <- stringr::str_match_all(posid, '([^\\+^\\-]+)([+-])(\\d+)')
  d <- data.frame(chr     = unlist(lapply(o, '[[', 2)),
                  strand  = unlist(lapply(o, '[[', 3)),
                  site    = as.integer(unlist(lapply(o, '[[', 4))),
                  subject = subject,
                  abund   = as.integer(abund))
  
  d$abundBins <- cut(d$abund, 10, labels = FALSE)
  d$site2 = d$site
  d$j = '0'
  
  file.create(file = outputFile)
  if(!is.na(position)) write(paste0('browser position ', position), file=outputFile, append = TRUE)
  
  write(paste0('track name="', title, ' positions" description="', title, '" colorByStrand="0,0,255 255,0,0"'), 
        file = outputFile, 
        append = TRUE)
  
  write.table(d[,c('chr', 'site', 'site2', 'subject', 'j', 'strand')], 
              file = outputFile, 
              sep = '\t', 
              append = TRUE, 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  
  write(paste0('track type="bedGraph" name="', title, 
               ' abund bins" color=0,128,0 visiblity=full autoScale=off viewLimits=0:10 maxHeightPixels=30:30:30'), 
        file=outputFile, 
        append = TRUE)
  
  write.table(d[,c('chr', 'site', 'site2', 'abundBins')], 
              file=outputFile, 
              sep = '\t', 
              append = TRUE, 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
}

