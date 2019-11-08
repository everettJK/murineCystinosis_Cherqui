createColorPalette <- function (n) 
{
  library(RColorBrewer)
  library(grDevices)
  colorRampPalette(brewer.pal(12, "Paired"))(n)
}


# UCSC visibility codes: 0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish. 
createUCSCintSiteTrack <- function(d, abundCuts = c(5,10,50), 
                                   posColors = c("#8C9DFF", "#6768E3", "#4234C7", "#1D00AB"),
                                   negColors = c("#FF8C8C", "#E35D5D", "#C72E2E", "#AB0000"),
                                   title='intSites', outputFile='track.ucsc', visibility = 1, 
                                   position = 'chr7:127471196-127495720', padSite = 0,
                                   siteLabel = NA){
  
  # Check function inputs.
  if(length(posColors) != length(negColors)) 
    stop('The pos and neg color vectors are not the same length.')

  if(length(abundCuts) != length(posColors) - 1) 
    stop('The number of aundance cut offs must be one less than the number of provided colors.')

  if(! all(c('start', 'end', 'strand', 'seqnames', 'estAbund') %in% names(d))) 
    stop("The expected column names 'start', 'end', 'strand', 'seqnames', 'estAbund' were not found.") 
  
  if(is.na(siteLabel) | ! siteLabel %in% names(d)) 
    stop('The siteLabel parameter is not define or can not be found in your data.')
  
  
  # Cut the abundance data. Abundance bins will be used to look up color codes.
  # We flank the provided cut break points with 0 and Inf in order to bin all values outside of breaks.
  cuts <- cut(d$estAbund, breaks = c(0, abundCuts, Inf), labels = FALSE)
  
  
  # Convert Hex color codes to RGB color codes. 
  # col2rgb returns a matrix, here we collapse the columns into comma delimited strings.
  posColors <- apply(grDevices::col2rgb(posColors), 2, paste0, collapse = ',')
  negColors <- apply(grDevices::col2rgb(negColors), 2, paste0, collapse = ',')
  
  
  # Create data fields needed for track table.
  d$score <- 0
  d$color <- ifelse(d$strand == '+', posColors[cuts], negColors[cuts])
  
  # Pad the site n NTs to increase visibility.
  if(padSite > 0){
    d$start <- floor(d$start - padSite/2)
    d$end   <- ceiling(d$end + padSite/2)
  }
  
  # Define track header.
  trackHead <- sprintf("track name='%s' description='%s' itemRgb='On' visibility=%s\nbrowser position %s",
                       title, title, visibility, position)
  
  # Write out track table.
  write(trackHead, file = outputFile, append = FALSE)
  write.table(d[, c('seqnames', 'start', 'end', siteLabel, 'score', 'strand', 'start', 'end', 'color')], 
              sep = '\t', col.names = FALSE, row.names = FALSE, file = outputFile, append = TRUE, quote = FALSE)
}
