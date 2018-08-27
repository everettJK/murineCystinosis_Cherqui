library(RMySQL)
group2.subjects <- scan(file = 'group2.subjects', what = 'character', sep = '\n')


# Read in sample data.
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="CYS"')


# Create a list of all GTSPs that passed through the INSPIIRED pipeline.
dbConn  <- dbConnect(MySQL(), group='intsites_miseq')
intSitesamples <- unname(unlist(dbGetQuery(dbConn, 'select sampleName from samples where sampleName like "%GTSP%"')))
intSitesamples <- unique(gsub('\\-\\d+$', '', intSitesamples))

# Only consider samples with intSite data.
samples <- samples[samples$SpecimenAccNum %in% intSitesamples,]

write(samples$Patient[! samples$Patient %in% group2.subjects], file = 'group1.subjects')