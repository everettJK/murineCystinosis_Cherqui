library(gt23)
library(RMySQL)
library(GenomicRanges)
library(dplyr)



stdIntSiteFragments <- function (frags, CPUs = 10, countsCol = "reads") 
{
  cluster <- parallel::makeCluster(CPUs)
  
  # parallel::clusterExport(cl = cluster, envir = environment(), 
  #                         varlist = c("countsCol"))
  
  frags <- gintools::standardize_sites(frags, counts.col = 'reads')
  
  
  frags <- unlist(GenomicRanges::GRangesList(parallel::parLapply(cluster, 
                                                                 as.list(split(frags, frags$sampleName)), function(x) {
                                                                   gintools::refine_breakpoints(x, counts.col = 'reads')
                                                                 })))
  parallel::stopCluster(cluster)
  frags$posid <- paste0(GenomicRanges::seqnames(frags), GenomicRanges::strand(frags), 
                        GenomicRanges::start(GenomicRanges::flank(frags, -1, 
                                                                  start = T)))
  ### dplyr::group_by(data.frame(frags), cellType, timePoint, start, 
  dplyr::group_by(data.frame(frags), sampleName, start, 
                  end, strand) %>% dplyr::mutate(reads = sum(reads)) %>% 
    dplyr::slice(1) %>% dplyr::ungroup() %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

getIntSiteData <-function (sampleDB.group, intSiteDB.group, patients = NULL, samples = NULL, roundTimePoints = TRUE, runIDs = NULL) 
{
  options(useFancyQuotes = FALSE)
  if (length(patients) > 0 && length(samples) > 0) 
    stop("Both patients and samples can not be defined.")
  if (length(patients) == 0 && length(samples) == 0) 
    stop("Either patients or samples must be defined.")
  dbConn1 <- DBI::dbConnect(RMySQL::MySQL(), group = sampleDB.group)
  dbConn2 <- DBI::dbConnect(RMySQL::MySQL(), group = intSiteDB.group)
  if (length(patients) > 0) {
    selectString <- "select SpecimenAccNum from gtsp where Patient='%s'"
    samples <- unname(unlist(sapply(patients, function(p) {
      unlist(DBI::dbGetQuery(dbConn1, sprintf(selectString, 
                                              p)))
    })))
  }
  if (length(samples) == 0) 
    stop("Error: no samples have been selected.")
  intSiteSamples <- DBI::dbGetQuery(dbConn2, "select * from samples")
  intSiteSamples$GTSP <- gsub("\\-\\d+$", "", intSiteSamples$sampleName)
  sampleIDs <- unique(base::subset(intSiteSamples, GTSP %in% 
                                     samples)$sampleID)
  if (length(sampleIDs) == 0) 
    stop("Error: no intSite sample ids have been selected.")
  replicateQuery <- paste("samples.sampleID in (", paste0(sampleIDs, 
                                                          collapse = ","), ")")
  if (!is.null(runIDs)) {
    runIDs <- paste0("(miseqid in (", paste0(sQuote(runIDs), 
                                             collapse = ","), ")) and ")
  }
  else {
    runIDs = ""
  }
  q <- sprintf("select position, chr, strand, breakpoint, count,\n               sampleName from sites left join samples on\n               sites.sampleID = samples.sampleID\n               left join pcrbreakpoints on\n               pcrbreakpoints.siteID = sites.siteID\n               where %s (%s)", 
               runIDs, replicateQuery)
  options(warn = -1)
  sampleData <- DBI::dbGetQuery(dbConn1, "select * from gtsp")
  sites <- dbGetQuery(dbConn2, q)
  options(warn = 0)
  if (nrow(sites) == 0) 
    return(GRanges())
  sites$GTSP <- as.character(sub("\\-\\d+", "", sites$sampleName))
  sites$patient <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 
                              "Patient"]
  sites$timePoint <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 
                                "Timepoint"]
  sites$cellType <- sampleData[match(sites$GTSP, sampleData$SpecimenAccNum), 
                               "CellType"]
  sites$timePoint <- toupper(sites$timePoint)
  sites$timePoint <- gsub("_", ".", sites$timePoint)
  sites$timePointType <- stringr::str_match(sites$timePoint, 
                                            "[DMY]")
  sites$timePointType[which(is.na(sites$timePointType))] <- "X"
  sites <- do.call(rbind, lapply(split(sites, sites$timePointType), 
                                 function(x) {
                                   n <- as.numeric(stringr::str_match(x$timePoint, "[\\d\\.]+"))
                                   if (x$timePointType[1] == "D") {
                                     x$timePointMonths <- n/30.4167
                                     x$timePointDays <- n
                                   }
                                   else if (x$timePointType[1] == "M") {
                                     x$timePointMonths <- n
                                     x$timePointDays <- n * 30.4167
                                   }
                                   else if (x$timePointType[1] == "Y") {
                                     x$timePointMonths <- n * 12
                                     x$timePointDays <- n * 365
                                   }
                                   else {
                                     message("Warning - could not determine date unit for: ", 
                                             x$timePointType[1])
                                     x$timePointMonths <- n
                                     x$timePointDays <- n
                                   }
                                   x
                                 }))
  sites$timePointType <- NULL
  if (roundTimePoints) 
    sites$timePointMonths <- base::round(sites$timePointMonths)
  GenomicRanges::GRanges(seqnames = S4Vectors::Rle(sites$chr), 
                         ranges = IRanges::IRanges(start = pmin(sites$position, 
                                                                sites$breakpoint), end = pmax(sites$position, sites$breakpoint)), 
                         strand = S4Vectors::Rle(sites$strand), reads = sites$count, 
                         patient = sites$patient, sampleName = sites$sampleName, 
                         GTSP = sites$GTSP, cellType = sites$cellType, timePoint = sites$timePoint, 
                         timePointDays = sites$timePointDays, timePointMonths = sites$timePointMonths)
}



calcReplicateAbundances <- function (gr) 
{
  if (!"posid" %in% names(mcols(gr))) {
    gr <- addPositionID(gr)
  }
  if (!"reads" %in% names(mcols(gr))) {
    stop("The provided GRange object does not have a 'reads' column")
  }
  if (!"GTSP" %in% names(mcols(gr))) {
    stop("The provided GRange object does not have a 'GTSP' column")
  }
  names(gr) <- NULL
  gr$startPosition <- start(flank(gr, -1, start = T))
  d <- as.data.frame(gr)
  d$frags <- 0
  d$w <- paste0(d$sampleName, "/", d$width)
  o <- dplyr::group_by(d, GTSP, posid) %>% dplyr::mutate(reads = sum(reads), 
                                                         frags = length(unique(w))) %>% dplyr::slice(1) %>% dplyr::ungroup() %>% 
    dplyr::mutate(estAbund = frags) %>% dplyr::group_by(GTSP) %>% 
    dplyr::mutate(relAbund = (estAbund/sum(estAbund)) * 100) %>% 
    dplyr::ungroup() %>% data.frame()
  o$start <- o$startPosition
  o$end <- o$startPosition
  o$w <- NULL
  o$startPosition <- NULL
  GenomicRanges::makeGRangesFromDataFrame(o, keep.extra.columns = TRUE)
}


frags1 <- getIntSiteData('specimen_management', 'intsites_miseq', samples = 'GTSP2737')
frags2 <- gt23::getDBgenomicFragments('GTSP2737', 'specimen_management', 'intsites_miseq')

length(frags1)
length(frags2)

stdFrags1 <- unlist(GRangesList(lapply(split(frags1, frags1$patient), function(p){
      gintools::standardize_sites(p, counts.col = 'reads')
    })))
length(stdFrags1)

stdFrags1 <- unlist(GRangesList(lapply(split(stdFrags1, stdFrags1$sampleName), function(p){
      gintools::refine_breakpoints(p, counts.col = 'reads')
    })))
length(stdFrags1)


stdFrags2 <- stdIntSiteFragments(frags2) 

length(stdFrags1)
length(stdFrags2)


sites1   <- calcReplicateAbundances(stdFrags1)
sites2   <- gt23::collapseReplicatesCalcAbunds(stdFrags2)

sum(sites1$estAbund)
sum(sites2$estAbund)

save.image('gt23.testing.RData')


