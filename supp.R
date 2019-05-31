



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

