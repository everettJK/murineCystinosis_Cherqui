options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(gt23)  # https://github.com/everettJK/package.geneTherapy.gt23
library(RMySQL)
library(parallel)
library(gtools)
library(GenomicRanges)
library(grDevices)
library(RColorBrewer)
library(tidyverse)
library(gintools)
source('./supp.R')
CPUs <- 40

savePointPrefix    <- 'group2'
reportSubjectsFile <- 'data/group2.subjects'
reportCellTransfersFile <- 'data/group2.cellTransfers_v2.tsv'


reportSubjects <- scan(reportSubjectsFile, what = 'character', sep = '\n')


# Read in supporting data files.
sample2organism <- read.table('data/sample2organism.tsv', sep='\t', strip.white = TRUE, header = TRUE)


# cellTransfers   <- read.table('data/cellTransfers.tsv', header=TRUE, sep='\t', check.names = FALSE)
cellTransfers <- data.frame(From = 'none', To = 'none')
if(file.exists(reportCellTransfersFile)) cellTransfers <- read.table(reportCellTransfersFile, header=TRUE, sep='\t', check.names = FALSE)


# Sanity check to make sure that each sample was processed against the correct genome.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
t(apply(sample2organism, 1, function(x){
  c(x, paste0(unique(dbGetQuery(dbConn, paste0('select * from samples where sampleName like "', x[1], '-%"'))$refGenome), collapse = ', '))
}))


# Read in sample data.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="CYS"')


# new
samples <- subset(samples, toupper(samples$Patient) %in% toupper(reportSubjects))

# Hot fix
### samples <- subset(samples, ! SpecimenAccNum %in% c('GTSP1708', 'GTSP1709'))


# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))


# Create a subject -> organism table
subjects <- do.call(rbind, lapply(split(samples, samples$Patient), function(x){
  data.frame(subject=x$Patient[1], organism=sample2organism[match(x$SpecimenAccNum, sample2organism$sample),]$organism)
}))


# Split subjects by organism, retrieve intSites and calculate intSite attributes.
intSites <- unlist(GRangesList(lapply(split(subjects, subjects$organism), function(x){
   # Retrieve intSites.
   intSites <- getIntSiteData('specimen_management', 'intsites_miseq', patients = x$subject)
   
   
   # Setup parallelization.
   cluster <- parallel::makeCluster(CPUs)
   
  
   intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$patient), function(p){
       gintools::standardize_sites(p, counts.col = 'reads')
   })))
   
    
   intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$sampleName), function(p){
      gintools::refine_breakpoints(p, counts.col = 'reads')
   })))
    
   
   # Merge replicate samples and estimate clonal abundances.
   intSites <- gt23::calcReplicateAbundances(intSites)
   
   # Create a splitting vector for parallelization.
   intSites$s <- ntile(seq_along(intSites), CPUs)
   
   # Add nearest gene and nearest oncogene annotations.
   names(intSites) <- NULL
   
   intSites$organism <- x$organism[1]
   intSites$genome   <- ifelse(x$organism[1] == 'human', 'hg38', 'mm9')
   
   intSites <- unlist(GRangesList(parLapply(cluster, split(intSites, intSites$s), function(p){
     source('./supp.R')
     library(dplyr)
     
     # Nearest gene boundary
     p$order <- 1:length(p)
     p <- nearestGenomicFeature(p, genome = p$genome[1])
     p <- p[order(p$order)]
     
     if(p$organism[1]== 'human'){
       geneList <- gt23::humanOncoGenesList
     } else {
       geneList <- gt23::mouseOncoGenesList
     }
     
     p$oncoGeneListLength <- length(geneList)
     
     # Nearest oncogene
     o <- nearestGenomicFeature(p, genome = p$genome[1], geneList = geneList)
     om <- data.frame(mcols(o))
     om <- om[order(om$order),]
     
     pm <- data.frame(mcols(p))
     pm <- pm[order(pm$order),]
     
     pm$nearestOncoFeature       <- om$nearestFeature
     pm$nearestOncoFeatureDist   <- om$nearestFeatureDist
     pm$nearestOncoFeatureStrand <- om$nearestFeatureStrand
     mcols(p) <- pm
    
     p
   })))
   
   stopCluster(cluster)
   
   intSites
})))


save.image(file = paste0('savePoints/', savePointPrefix, '.1.RData'))

samples[which(samples$Timepoint == 0),]$Timepoint <- "D0"
samples$Timepoint <- toupper(samples$Timepoint)
samples[which(samples$Timepoint == '6M'),]$Timepoint <- "M6"
intSites$timePoint <- toupper(intSites$timePoint)
intSites$timePoint <- gsub('6M', 'M6', intSites$timePoint)
intSites$timePoint <- gsub('^0$', 'D0', intSites$timePoint)

intSites$cellType <- gsub('BM\\s+GM', 'BM Myeloid cells', intSites$cellType)
intSites$cellType <- gsub('BM\\s+B\\s+Cells', 'B-cells', intSites$cellType)
intSites$cellType <- gsub('BM\\s+B\\-Cells', 'B-cells', intSites$cellType)


d <- data.frame(mcols(intSites))
d[which(d$GTSP == 'GTSP2347'),]$cellType <- 'Hematoma'
mcols(intSites) <- d



# Hot fixes
r <- unique(scan(file = 'data/group2_required_samples', what = 'character'))

#  r[! r %in% intSitesamples]
# "GTSP2348" "GTSP2327"


# Add VCN values.
intSites$VCN <- sapply(intSites$GTSP, function(x){ round(samples[match(x, samples$SpecimenAccNum),]$VCN, digits=3) })
if(length(which(intSites$VCN == 0) > 0)) intSites[which(intSites$VCN == 0)]$VCN <- NA


# First, check with CYS samples are in the full list of INSPIIRED samples and then determine which 
# of those samples are not in the intSite object which requires at least 1 site to be found for inclussion. 
processedSamples <- samples$SpecimenAccNum[samples$SpecimenAccNum %in% intSitesamples]
samplesNoIntSitesFound <- processedSamples[!processedSamples %in% intSites$GTSP]

failedSampleTable <-
  samples %>%
  select(SpecimenAccNum, CellType, Patient, Timepoint, SpecimenInfo) %>%
  filter(SpecimenAccNum %in% samplesNoIntSitesFound) %>%
  mutate(SpecimenInfo = ifelse(SpecimenInfo == 'Mouse', 'none', SpecimenInfo))


# Create an organism specific effort table.
summaryTable <- 
  intSites %>%
  data.frame() %>%
  mutate(patientPosid = paste(patient, posid)) %>%
  group_by(organism) %>%
  summarise(samples = n_distinct(GTSP),
            nReads  = ppNum(sum(reads)),
            nFrags  = ppNum(sum(estAbund)),
            nSites  = ppNum(n_distinct(patientPosid))) %>%
  ungroup()


# Create a read depth visualization.
intSiteReadsPlot <-
  intSites %>%
  data.frame() %>%
  group_by(GTSP, posid) %>%
  summarise(nReads = sum(reads),
            group  = ifelse(organism == 'human', 'Human', 
                            ### ifelse(patient %in% cellTransfers$From, 'Mouse donor',
                            ifelse(GTSP %in% cellTransfers$From, 'Mouse donor',
                                   ### ifelse(patient %in% cellTransfers$To, 'Mouse recipient', 'Mouse')))) %>%
                                   ifelse(GTSP %in% cellTransfers$To, 'Mouse recipient', 'Mouse')))) %>%
  ungroup() %>%
  arrange(group) %>%
  mutate(GTSP = factor(GTSP, levels = unique(GTSP))) %>%
  ggplot(aes(GTSP, log10(nReads), fill=group)) + 
    theme_bw() +
    geom_point(shape=22, alpha=0.05, stroke = 0, size=4) + 
    coord_flip() +
    guides(fill = guide_legend(title="Sample type", override.aes = list(alpha = 1))) +
    scale_fill_manual(values = c('red', 'green', 'blue', 'gold3')) +
    labs(y = 'log10(number of reads)', x = 'Sample')


# Create a table of human subjects with the percent of sites near suspect oncogenes.
# humanSitesNearOnco <- 
#   data.frame(subset(intSites, organism=='human')) %>%
#   group_by(patient) %>%
#   summarise(percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
#   ungroup() %>%
#   arrange(percentNearOnco) %>%
#   mutate(source = 'CYS') %>%
#   data.frame()


# Read in previously published WAS trial d0 data and rename the subjects and then determin the percentage
# of intSites near oncogenes.

# n <- 1
# WASintSites_d0 <- readRDS('data/WAS_d0_intSites.rds')
# WASintSites_d0 <- unlist(GRangesList(lapply(split(WASintSites_d0, WASintSites_d0$patient), 
#                                             function(x){ x$patient <- paste('WAS subject', n); n <<- n+1; x})))
# 
# 
# wasSitesNearOnco <- 
#   data.frame(WASintSites_d0) %>%
#   group_by(patient) %>%
#   summarise(percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
#   ungroup() %>%
#   arrange(percentNearOnco) %>%
#   mutate(source = 'WAS') %>%
#   data.frame()
# 
# 
# # Create a bar plot comparing the percentage of intSites near oncogenes for human subjects vs WAS d0 subjects.
# WASvsHumanSubjects <- 
#   bind_rows(humanSitesNearOnco, wasSitesNearOnco) %>%
#   arrange(percentNearOnco) %>%
#   mutate(patient = factor(patient, levels=unique(patient))) %>%
#   ggplot(aes(patient, percentNearOnco, fill=source)) +
#     theme_bw() +
#     scale_fill_manual(values=c('gray75', 'blue')) +
#     geom_bar(stat='identity') +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     guides(fill = FALSE) +
#     scale_y_continuous(limits = c(0, 0.33), labels = scales::percent) +
#     labs(x='Subject', y='Sites near oncogenes')
    

# Create similiar data frames and plots for the mouse subject by comparing the mouse subjects to a previously 
# published mouse study.

mouseSitesNearOnco <- 
  data.frame(subset(intSites, organism=='mouse')) %>%
  ### filter(patient %in% cellTransfers$From) %>%
  filter(GTSP %in% cellTransfers$From) %>%
  group_by(patient) %>%
  summarise(percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
  ungroup() %>%
  filter(percentNearOnco > 0) %>%
  arrange(percentNearOnco) %>%
  mutate(source = 'CYS') %>%
  data.frame()


load('data/PMC3129560_mouseTrial/PMC3129560.RData')
PMC3129560.intSiteData$patient <- 'PMC3129560'
PMC3129560.intSiteData <- gt23::addPositionID(PMC3129560.intSiteData)


PMC3129560SitesNearOnco <- 
  data.frame(PMC3129560.intSiteData) %>%
  summarise(patient=patient[1],
            percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
  filter(percentNearOnco > 0) %>%
  mutate(source = 'PMC3129560') %>%
  data.frame()


PMC3129560vsMouseSubjects <-
  bind_rows(mouseSitesNearOnco, PMC3129560SitesNearOnco) %>%
  arrange(percentNearOnco) %>%
  mutate(patient = factor(patient, levels=unique(patient))) %>%
  ggplot(aes(patient, percentNearOnco, fill=source)) +
  theme_bw() +
  scale_fill_manual(values=c('gray75', 'blue')) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE) +
  scale_y_continuous(limits = c(0, 0.10), labels = scales::percent) +
  labs(x='Subject', y='Sites near oncogenes')


# Create a series of distributions showing the integrations sites vs the number of 
# associated fragments (inferred cells) for the study samples as well as the previously 
# published WAS time points.

# prevWASd0IntSiteFrags <- 
#   WASintSites_d0 %>%
#   data.frame() %>%
#   group_by(estAbund, timePoint) %>%
#   summarise(source = 'WAS', nSites = n_distinct(posid)) %>%
#   ungroup() 
# 
# prevWASintSiteFrags <- 
#   readRDS('data/prevWASintSites.rds') %>%
#   data.frame() %>%
#   group_by(estAbund, timePoint) %>%
#   summarise(source = 'WAS', nSites = n_distinct(posid)) %>%
#   ungroup() 
# 
# cysIntSiteFrags <- 
#   intSites %>%
#   data.frame() %>%
#   filter(organism == 'human') %>%
#   group_by(estAbund, timePoint) %>%
#   summarise(source = 'CYS', nSites = n_distinct(posid)) %>%
#   ungroup() 
# 
# intSiteFragPlot <-
#   rbind(prevWASd0IntSiteFrags, prevWASintSiteFrags, cysIntSiteFrags) %>% 
#   filter(timePoint %in% c('d0', 'D14', 'm6', 'm12')) %>%
#   mutate('Data set' = toupper(paste(source, timePoint))) %>%
#   ggplot(aes(estAbund, log2(nSites+1), fill=`Data set`)) +
#     theme_bw() +
#     geom_bar(stat='identity') +
#     scale_fill_manual(values=c('green3', 'dodgerblue', 'gold1', 'red')) +
#     xlim(c(1, 100)) +
#     labs(x='Clones', y='log2(Number of integration sites + 1)')


# Create intSite heat maps.

# # Create a list of chromosome lengths for select chromosomes, ie. [['chr1']] <- 248956422.
# library(BSgenome.Hsapiens.UCSC.hg38)
# names(intSites)   <- NULL
# names(WASintSites_d0) <- NULL
# chromosomeLengths <- sapply(rev(paste0("chr", c(seq(1:21), "X", "Y"))),
#                             function(x){length(BSgenome.Hsapiens.UCSC.hg38[[x]])},
#                             simplify = FALSE, USE.NAMES = TRUE)
# 
# humanIntSiteMap <- intSiteDistributionPlot(subset(intSites, organism == 'human'), chromosomeLengths, alpha = 0.025)
# WASintSites_d0_map <- intSiteDistributionPlot(WASintSites_d0, chromosomeLengths, alpha = 0.2)


library(BSgenome.Mmusculus.UCSC.mm9)
chromosomeLengths <- sapply(rev(paste0("chr", c(seq(1:19), "X", "Y"))),
                            function(x){length(BSgenome.Mmusculus.UCSC.mm9[[x]])},
                            simplify = FALSE, USE.NAMES = TRUE)

names(intSites) <- NULL
mouseIntSiteMap <- intSiteDistributionPlot(subset(intSites, organism == 'mouse'), chromosomeLengths, alpha = 0.2)
### mouseDonorIntSiteMap <- intSiteDistributionPlot(subset(intSites, patient %in% cellTransfers$From), chromosomeLengths, alpha = 0.3)
mouseDonorIntSiteMap <- intSiteDistributionPlot(subset(intSites, GTSP %in% cellTransfers$From), chromosomeLengths, alpha = 0.3)
### mouseRecipientIntSiteMap <- intSiteDistributionPlot(subset(intSites, patient %in% cellTransfers$To), chromosomeLengths, alpha = 0.5)
mouseRecipientIntSiteMap <- intSiteDistributionPlot(subset(intSites, GTSP %in% cellTransfers$To), chromosomeLengths, alpha = 0.5)


# Create relative abundance plots for the human samples and store them as a list of grobs so that they 
# can be arranged in the report.

o <-
  intSites %>%
  data.frame() %>%
  filter(organism == 'mouse') %>%
  group_by(GTSP) %>%
  mutate(nSites = n_distinct(posid)) %>%
  arrange(desc(relAbund)) %>%
  filter(between(row_number(), 1, 25)) %>%
  select(GTSP, patient, timePoint, cellType, relAbund, nearestFeature, posid, nSites) %>%
  do(add_row(.,
             GTSP=.$GTSP[1],
              patient=.$patient[1],
              relAbund=100-sum(.$relAbund),
              nearestFeature='LowAbund',
              posid='x',
              nSites=.$nSites[1],
              .before = 1)) %>%
  ungroup()

mouseRelAbundPlots <-
  lapply(split(o, o$GTSP), function(x){
    x$nearestFeature <- factor(x$nearestFeature, levels = unique(x$nearestFeature))
    ggplot(x, aes(GTSP, relAbund/100, fill = nearestFeature)) +
      theme_bw() +
      geom_bar(stat='identity') +
      scale_fill_manual(values = c('gray90', colorRampPalette(brewer.pal(12, "Paired"))(25))) +
      labs(x='', y='') +
      scale_y_continuous(labels = scales::percent) +
      theme(legend.position="none") +
      ggtitle(paste0(x$patient[1], '\n', x$cellType[nrow(x)], '\n', x$timePoint[nrow(x)], ', ', ppNum(x$nSites[1]), ' sites')) +
      theme(plot.title = element_text(size = 6.5)) +
      theme(plot.margin = unit(c(0,0,0,0), "cm"))
  })

rm(o)


# Create a table of clones that exceeded 20% relative abundance.
abundantClones20 <-
  intSites %>%
  data.frame() %>%
  select(patient, organism, timePoint, cellType, posid, relAbund, estAbund, nearestFeature) %>%
  filter(relAbund > 20) %>%
  arrange(organism) %>%
  mutate(relAbund = sprintf("%.02f%%", relAbund))


# Create a sample summary of all analyzed samples.
sampleSummary <-
  intSites %>%
  data.frame() %>%
  select(patient, GTSP, organism, timePoint, cellType, posid, relAbund, estAbund, nearestFeature, VCN) %>%
  group_by(organism, GTSP) %>%
  summarise(Subject = patient[1],
            'Cell type' = cellType[1],
            'VCN' = VCN[1],
            'Time point' = timePoint[1],
            'Number inferred cells' = ppNum(sum(estAbund)),
            'Number of intSites' = ppNum(length(unique(posid)))) %>%
  ungroup()
  
  
# Currently the returned intSite object contains 'TRUE' or NA -- NA breaks ifelse.
intSites$inFeature[is.na(intSites$inFeature)] <- 'FALSE'

intSites$nearestFeature2 <- 
  intSites %>%
  data.frame() %>%
  mutate(nearestFeature2 = paste0(nearestFeature, ' ')) %>% 
  mutate(nearestFeature2 = ifelse(inFeature == 'TRUE', paste0(nearestFeature2, '*'), nearestFeature2)) %>%
  mutate(nearestFeature2 = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(nearestFeature2, '~'), nearestFeature2)) %>%
  select(nearestFeature2) %>%
  unlist() %>%
  unname()


# Cycle through the cell transplant trials and create relative abunance plots, matrix of intSites
# near oncogenes, and data frames of intSites that persist in the recipient mice.


emptyRecipientPlotLabels <- list('pCN774' = 'pCN948\n(no sites)', 'pCN809' = 'pCN952\n(no sites)')


transferTrials <- lapply(1:nrow(cellTransfers), function(i){
  d <- cellTransfers[i,]
  
  ## a <- data.frame(subset(intSites, patient == d$From))
  ## b <- data.frame(subset(intSites, patient == d$To))
  
  #if(d$From == 'GTSP1697') browser()
  
  a <- data.frame(subset(intSites, GTSP == d$From))
  b <- data.frame(subset(intSites, GTSP == d$To))
  
  createPlotData <- function(x){
    arrange(x, desc(relAbund)) %>%
    mutate(label1 = paste0(x$patient[1], '\n',
                           x$cellType[1], '\n',
                          'VCN: ', VCN[1], '\n',
                          'Unique sites: ', ppNum(length(unique(posid))), '\n',
                          'Inferred cells: ', numShortHand(sum(estAbund)), '\n',
                          x$timePoint[1], ' / ', x$cellType[1])) %>%
    mutate(label2 = paste0(nearestFeature2, '\n', posid)) %>%
    filter(between(row_number(), 1, 12)) %>%
    select(label1, label2, relAbund) %>%
    add_row(relAbund=100-sum(.$relAbund),
            label1 = .$label1[1],
            label2 = 'LowAbund',
            .before = 1)
  }
  
  plotData <- bind_rows(createPlotData(a), createPlotData(b))
  
  if(a$patient[1] %in% names(emptyRecipientPlotLabels) & length(which(is.na(plotData$label1))) > 0){
    plotData[which(is.na(plotData$label1)),]$label1 <- emptyRecipientPlotLabels[[a$patient[1]]]
  }
  
  plot <- 
    plotData %>%
    mutate(label1 = factor(label1, levels=unique(label1))) %>%
    mutate(label2 = factor(label2, levels = unique(label2))) %>%
    mutate(label2 = fct_relevel(label2, 'LowAbund')) %>%
    arrange(desc(relAbund)) %>%
    ggplot(aes(label1, relAbund, fill=label2)) +
      theme_bw() +
      geom_bar(stat='identity') +
      scale_fill_manual(name = 'intSites', values = c('gray90', createColorPalette(24))) +
      labs(x='', y='Relative abundance') +
    guides(fill=guide_legend(ncol=2)) +
    theme(legend.key.size = unit(2, "line"), legend.text=element_text(size=10))
  
  guides(shape = guide_legend(override.aes = list(size = 5)))
  
  
  
  m <- matrix(c(sum(abs(a$nearestOncoFeatureDist) > 50000, na.rm = TRUE),  sum(abs(a$nearestOncoFeatureDist) <= 50000, na.rm = TRUE),
                sum(abs(b$nearestOncoFeatureDist) > 50000, na.rm = TRUE),  sum(abs(b$nearestOncoFeatureDist) <= 50000, na.rm = TRUE)), 
              byrow = TRUE, 
              nrow = 2,
              dimnames = list(c(a$patient[1], b$patient[1]), c('Not near onco', 'Near onco')))

  
  sharedSites <- bind_rows(lapply(b$posid[unique(b$posid) %in% unique(a$posid)], function(posID){
    data.frame(Donor = a$patient[1],
               Recipient = b$patient[1],
               intSite = posID,
               'Donor cells' = ppNum(sum(subset(a, posid == posID)$estAbund)),
               'Recipient cells' = ppNum(sum(subset(b, posid == posID)$estAbund)),
               check.names = FALSE) }))
  
  list(plot = plot, m = m, sharedSites = sharedSites)
})


# Assemble the intSite persistence table
cellTransfer_intSites_table <- bind_rows(lapply(transferTrials, '[[', 3))


# Create UCSC track files.
# o <- subset(intSites, organism == 'human')
# createUCSCintSiteAbundTrack(o$posid, o$estAbund, subject = 'CYS_human', title = 'CYS_human', outputFile = 'UCSC_CYS_human.group1.ucsc')
# system(paste0('scp UCSC_CYS_human.group1.ucsc  microb120:/usr/share/nginx/html/UCSC/cherqui/'))
# file.remove('UCSC_CYS_human.group1.ucsc')


o <- subset(intSites, organism == 'mouse')
createUCSCintSiteAbundTrack(o$posid, o$estAbund, subject = 'CYS_mouse', title = 'CYS_mouse', outputFile = 'UCSC_CYS_mouse.group2.ucsc')
system(paste0('scp UCSC_CYS_mouse.group2.ucsc  microb120:/usr/share/nginx/html/UCSC/cherqui/'))
file.remove('UCSC_CYS_mouse.group2.ucsc')


# Report shortcuts.
#humanGenomePercentOnco <- round((length(gt23::humanOncoGenesList)) / length(unique(gt23::hg38.refSeqGenesGRanges$name2))*100, digits=2)
mouseGenomePercentOnco <- round((length(gt23::mouseOncoGenesList)) / length(unique(gt23::mm9.refSeqGenesGRanges))*100, digits=2)


# Save data for report generation.
save.image(file='project.group2.RData')


# # Patient check
# p <- scan('group2.check', what = 'character', sep = '\n')
# i <- unique(c(intSites$GTSP, failedSampleTable$SpecimenAccNum))
# s <- unique(c(intSites$patient, failedSampleTable$Patient))
# 
# # Are all the patients in the patient check list accounted for in the data?
# table(p %in% i)
# 
# # Are there any patients in the data not in the check list?
# table(i %in% p)


