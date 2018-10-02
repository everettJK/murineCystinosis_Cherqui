d <- read.table('VCN_value_updates', sep = '\t', header = TRUE)


library(RMySQL)
invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')

apply(d, 1, function(x){
  #browser()
  q <- paste0('update gtsp set VCN = "', x[2], '" where SpecimenAccNum = "', x[1], '"')
  message(q)
  dbSendQuery(dbConn, q)
  q <- paste0('select VCN from  gtsp where SpecimenAccNum = "', x[1], '"')
  r <- dbGetQuery(dbConn, q)
  message(as.character(r), ' == ', x[2])
})
