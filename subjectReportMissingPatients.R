options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
library(RMySQL)

p <- unique(c(scan(what = 'character', file = 'data/group1.subjects'), scan(what = 'character', file = 'data/group2.subjects')))
r <- read.table(file = '../project.geneTherapy.subjectReport/subjectParameters.full', sep = ';', header = TRUE)

# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))


# Read in sample data.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="CYS"')

samples <- subset(samples, SpecimenAccNum %in% intSitesamples)

missingPatients <- p[! p %in% r$Patient]

write(paste0('"CYS;"', missingPatients[missingPatients %in% samples$Patient], '";"";"mm9"'), file = 'subjectReportPatients')



