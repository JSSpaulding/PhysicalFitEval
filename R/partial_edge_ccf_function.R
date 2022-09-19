partial_edge_ccf <- function(edge1, edge2, pb=TRUE){
  ## Align at Xmin ----
  if(edge1$x[1] < edge2$x[1]){
    shift <- edge2$x[1] - edge1$x[1]
    edge2$x <- edge2$x - shift
  }
  if(edge1$x[1] > edge2$x[1]){
    shift <- edge1$x[1] - edge2$x[1]
    edge1$x <- edge1$x - shift
  }
  ## Determine Lag/Movement of Partial Along Edge 1 ----
  lags <- seq(0,max(edge1$x)-max(edge2$x,na.rm = T))
  ## CCF of Partial Along Edge 1 ----
  jj <- 1
  corresult <- c()
  if(pb==TRUE){
    pb = txtProgressBar(min = 0, max = max(lags), initial = 0, style = 3) 
    for (i in lags) {
      setTxtProgressBar(pb,i)
      #Shift Edge 2 across Edge 1 --
      e2 <- edge2
      e2$x <- edge2$x + i
      #Determine Comparison Area
      comp_area <- c(min(e2$x),max(e2$x))
      e1 <- edge1[edge1$x > comp_area[1] & edge1$x < comp_area[2],]
      #Pad Edge 1 --
      e1int <- suppressWarnings(approx(e1$x,e1$y,xout = e2$x))
      e1int <- rbind(e1,e1int)
      edge1a <- data.frame(x=e1int$x,
                           y=e1int$y)
      edge1a <- edge1a %>% 
        # Order data by x
        arrange(x) %>%
        # Subset the coord columns
        select(x, y)
      #Pad Edge 2 --
      e2int <- suppressWarnings(approx(e2$x,e2$y,xout = e1$x))
      e2int <- rbind(e2,e2int)
      edge2a <- data.frame(x=e2int$x,
                           y=e2int$y)
      edge2a <- edge2a %>% 
        # Order data by x
        arrange(x) %>%
        # Subset the coord columns
        select(x, y)
      cors <- cor.test(edge1a$y,edge2a$y)
      corresult[jj] <- cors$estimate
      jj <- jj + 1
    }
  } else{
    for (i in lags) {
      #Shift Edge 2 across Edge 1 --
      e2 <- edge2
      e2$x <- edge2$x + i
      #Determine Comparison Area
      comp_area <- c(min(e2$x),max(e2$x))
      e1 <- edge1[edge1$x > comp_area[1] & edge1$x < comp_area[2],]
      #Pad Edge 1 --
      e1int <- suppressWarnings(approx(e1$x,e1$y,xout = e2$x))
      e1int <- rbind(e1,e1int)
      edge1a <- data.frame(x=e1int$x,
                           y=e1int$y)
      edge1a <- edge1a %>% 
        # Order data by x
        arrange(x) %>%
        # Subset the coord columns
        select(x, y)
      #Pad Edge 2 --
      e2int <- suppressWarnings(approx(e2$x,e2$y,xout = e1$x))
      e2int <- rbind(e2,e2int)
      edge2a <- data.frame(x=e2int$x,
                           y=e2int$y)
      edge2a <- edge2a %>% 
        # Order data by x
        arrange(x) %>%
        # Subset the coord columns
        select(x, y)
      cors <- cor.test(edge1a$y,edge2a$y)
      corresult[jj] <- cors$estimate
      jj <- jj + 1
    }
  }
  return(max(abs(corresult)))
}
