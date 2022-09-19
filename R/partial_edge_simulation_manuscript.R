#### Characterization of Fracture Match Associations with Auto. Image Processing
#### Supplemental Materials
#### Physical Fit Comparisons

#### REQUIRED PACKAGES & DIRECTORIES ===========================================
maindir <- "D:/Workspace/PhysicalFitEval"
duct_samples <- paste0(maindir, "/tape-samples")
#outdir <- paste0(maindir, "/results")

setwd(maindir)
source("partial_edge_ccf_function.R") #Load partial edge comparison script as function

## Packages --
pkgs <- c("imager","dplyr","pbapply","stringr")
install.packages(pkgs) #installs packages listed in pkgs
lapply(pkgs, library, character.only = TRUE) #loads required packages in pkgs


#### LOAD FILES & GENERATE PARTIAL EDGES =======================================
setwd(duct_samples)
files <- list.files() #lists all file names in directory
q_edges <- seq(1, length(list.files()),2) #define questioned edges
k_edges <- seq(2, length(list.files()),2) #define known edges

partial_edges <- list()
p_edge_names <- c()
pb = txtProgressBar(min = 0, max = length(q_edges), initial = 0, style = 3) 
for (i in 1:length(q_edges)) {
  setTxtProgressBar(pb,i)
  ## STEP 1 - Read in Questioned Edge Image
  q_name <- files[q_edges[i]]
  samp1 <- load.image(q_name)
  #Rotate 90 degrees
  samp1 <- imrotate(samp1, 90)
  samp1 <- threshold(samp1, thr = "auto", approx = TRUE, adjust = 1.5)
  #Edge Detection
  edge1 <- tryCatch(cannyEdges(samp1))
  out1 <- data.frame(which(edge1, arr.ind = TRUE))
  names(out1) <- c('x','y','a','b')
  out1 <- out1 %>% 
    #Order data by x
    arrange(x) %>%
    #Subset the coord columns
    select(x, y)
  
  ## STEP 2 - Read in Known Edge Image
  k_name <- files[k_edges[i]]
  samp2 <- load.image(k_name)
  #Rotate 90 degrees
  samp2 <- imrotate(samp2, 90)
  samp2 <- threshold(samp2, thr = "auto", approx = TRUE, adjust = 1.5)
  #Edge Detection
  edge2 <- cannyEdges(samp2)
  out2 <- data.frame(which(edge2, arr.ind = TRUE))
  names(out2) <- c('x','y','a','b')
  out2 <- out2 %>% 
    #Order data by x
    arrange(x) %>%
    #Subset the coord columns
    select(x, y)
  
  ## STEP 3 - Align at Xmin
  if(out1$x[1] < out2$x[1]){
    shift <- out2$x[1] - out1$x[1]
    out1_align <- out1
    out2_align <- out2
    out2_align$x <- out2$x - shift
  } else{
    shift <- out1$x[1] - out2$x[1]
    out1_align <- out1
    out2_align <- out2
    out1_align$x <- out1$x - shift
  }
  
  ## STEP 4 - Translation of Known Edge: Rotation
  out1_crop <- out1_align[out1_align$x > min(out1_align$x) + 
                            100 & out1_align$x < max(out1_align$x) - 100 ,]
  out2_crop <- out2_align[out2_align$x > min(out2_align$x) + 
                            100 & out2_align$x < max(out2_align$x) - 100 ,]
  #Create Linear Model of Edge for Rotation
  out1_lm <- lm(y ~ x, data=out1_crop)
  out2_lm <- lm(y ~ x, data=out2_crop)
  out1_a <- out1_lm$coefficients[2]
  out2_a <- out2_lm$coefficients[2]
  out2_align <- data.frame(spdep::Rotation(out2_align, out1_a - out2_a))
  names(out2_align) <- c('x','y')
  
  ## STEP 5 - Round to Nearest Integer Such as Scan Data
  out2_align$x <- round(out2_align$x)
  out2_align$y <- round(out2_align$y)
  
  ## STEP 6 - Incremental Overlay of Edges: Find optimal alignment on Y-axis
  out2_align$y <- out2_align$y - 150
  overlap <- NULL
  jj <- 1
  for (j in 1:300) {
    shifted <- data.frame(x=out2_align$x,y=out2_align$y + j)
    overlap[jj] <- nrow(intersect(out1_align,shifted))
    jj <- jj + 1
  }
  out2_align$y <- out2_align$y + which.max(overlap)
  
  ## STEP 7 - Crop Various Sizes of Each Edge
  #Edge 1 of Mated Pair
  p_edges1 <- list()
  jj <- 1
  for (k in seq(min(out1_align$x),max(out1_align$x),200)) {
    p_edges1[[jj]] <- out1_align[out1_align$x <= k ,]
    jj <- jj + 1
  }
  p_edges1[[jj]] <- out1_align
  if(length(partial_edges)==0){
    partial_edges[[1]] <- p_edges1
  } else{
    partial_edges <- rlist::list.append(partial_edges,p_edges1)
  }
  #Name of edge1
  if(length(p_edge_names)==0){
    p_edge_names[1] <- gsub("\\..*","",q_name)
  } else{
    p_edge_names <- rlist::list.append(p_edge_names,gsub("\\..*","",q_name))
  }
  #Edge 2 of Mated Pair
  p_edges2 <- list()
  jj <- 1
  for (k in seq(min(out2_align$x),max(out2_align$x),200)) {
    p_edges2[[jj]] <- out2_align[out2_align$x <= k ,]
    jj <- jj + 1
  }
  p_edges2[[jj]] <- out2_align
  partial_edges <- rlist::list.append(partial_edges,p_edges1)
  
  #Name of edge2
  p_edge_names <- rlist::list.append(p_edge_names,gsub("\\..*","",k_name))
}

#Rename partial edges (organization)
names(partial_edges) <- p_edge_names

#### DEFINE PAIRWISE COMPARISONS ===============================================
comparisons <- expand.grid(p_edge_names,files)
#Use regular expressions to find self comparisons
matches <- c()
for (i in 1:nrow(comparisons)) {
  a <- comparisons[i,1]
  b <- gsub("\\..*","",comparisons[i,2])
  if(a==b){
    if(length(matches)==0){
      matches[1] <- i
    } else{
      matches <- rlist::list.append(matches,i)
    }
  }
}
comparisons <- comparisons[-c(matches),] #remove self comparisons

#### PARTIAL EDGE COMPARISON LOOP ==============================================
# Note: This loop, with all samples is computationally demanding with a run time
# of several weeks. To complete, this loop was run in a computer lab and with
# segmentation for manual parallelization to facilitate computation.

comb_out <- list()
pb = txtProgressBar(min = 0, max = nrow(comparisons), initial = 0, style = 3) 
for (i in 1:nrow(comparisons)) {
  edges_1 <- partial_edges[[comparisons[i,1]]]
  edge1_full <- edges_1[[20]]
  edges_2 <- partial_edges[[comparisons[i,2]]]
  edge2_full <- edges_2[[20]]
  
  for (j in 1:length(edges_2)) {
    #Compare all E1 with E2; Return width %, CCF, lag position (for RMP)
    edge1 <- edge1_full
    edge2 <- edges_2[[j]]
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
    # widthpercent <- c()
    # ccf <- c()
    # maxpos <- c()
    #pb = txtProgressBar(min = 0, max = length(lags), initial = 0, style = 3) 
    for (k in lags) {
      #setTxtProgressBar(pb,jj)
      #Shift Edge 2 across Edge 1 --
      e2 <- edge2
      e2$x <- edge2$x + k
      #Determine Comparison Area
      comp_area <- c(min(e2$x),max(e2$x))
      e1 <- edge1[edge1$x > comp_area[1] & edge1$x < comp_area[2],]
      
      if(nrow(e1)==0){
        corresult[jj] <- NA
      } else{
        #Pad Edge 1 --
        #if only 1 x value, simulate another b/c 2 needed
        if(length(unique(e1$x))==1){
          e1[nrow(e1) + 1,] = list(unique(e1$x)+1,e1$y[nrow(e1)]) 
        }
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
        #if only 1 x value, simulate another b/c 2 needed
        if(length(unique(e2$x))==1){
          e2[nrow(e2) + 1,] = list(unique(e2$x)+1,e2$y[nrow(e2)])
        }
        e2int <- suppressWarnings(approx(e2$x,e2$y,xout = e1$x))
        e2int <- rbind(e2,e2int)
        edge2a <- data.frame(x=e2int$x,
                             y=e2int$y)
        edge2a <- edge2a %>% 
          # Order data by x
          arrange(x) %>%
          # Subset the coord columns
          select(x, y)
        cors <- suppressWarnings(cor.test(edge1a$y,edge2a$y))
        corresult[jj] <- cors$estimate
      }
      jj <- jj + 1
    }
    if(length(comb_out)==0){
      comb_out <- data.frame(full_section=comparisons[i,1],
                             partial_section=comparisons[i,2],
                             width=round(nrow(edges_1[[j]]) / nrow(edge1_full), 
                                         digits = 2),
                             ccf=max(abs(corresult)),
                             position=lags[which.max(abs(corresult))])
    } else{
      temp_dat <- data.frame(full_section=comparisons[i,1],
                             partial_section=comparisons[i,2],
                             width=round(nrow(edges_1[[j]]) / nrow(edge1_full),
                                         digits = 2),
                             ccf=max(abs(corresult)),
                             position=lags[which.max(abs(corresult))])
      comb_out <- rbind(comb_out,temp_dat)
    }
    
  }
  
}


#### SORT RESULTS INTO KNOWN MATCHES AND KNOWN NON-MATCHES =====================
pb = txtProgressBar(min = 0, max = nrow(comb_out), initial = 0, style = 3) 
for (i in 1:nrow(comb_out)) {
  setTxtProgressBar(pb,i)
  if(str_sub(comb_out$full[i],1,-2) == str_sub(comb_out$partial[i],1,-2)){
    comb_out$type[i] <- "KM"
  } else{
    comb_out$type[i] <- "KNM"
  }
}
comb_out$ccf <- abs(comb_out$ccf)

KM <- subset(comb_out,comb_out$type=="KM")
KM_s <- KM[order(KM$width),]
KNM <- subset(comb_out,comb_out$type=="KNM")
KNM_s <- KNM[order(KNM$width),]


#### COLLATE RESULTS BY EDGE WIDTH =============================================
avg <- NULL
percent <- unique(KM_s$width)
jj <- 1
for (i in percent) {
  widthset <- subset(KM_s, KM_s$width==i)
  avg[jj] <- mean(widthset$ccf, na.rm = TRUE)
  jj <- jj + 1
}
comb <- cbind(percent,avg)
comb <- comb[order(percent),]


#### KM & KNM CCF BY EDGE WIDTH ================================================
#Figure 7 construction
plot(comb[,1], comb[,2], main = "Combined Duct Edge CCF Plots", 
     lwd = 2, ylim = c(0,1), xlab = "Edge Percent (%)", ylab = "CCF",pch=16)
#Trend for KM
comb <- data.frame(comb)
plx <- predict(loess(comb$avg ~ comb$percent, span = 0.3), se = T)
l <- msir::loess.sd(comb$avg ~ comb$percent, nsigma = 1.96)

#Trend for KNM
KNM_avg <- NULL
KNM_percent <- unique(KNM_s$width)
jj <- 1
for (i in KNM_percent) {
  widthset <- subset(KNM_s, KNM_s$width==i)
  KNM_avg[jj] <- mean(widthset$ccf, na.rm = TRUE)
  jj <- jj + 1
}
KNM_comb <- cbind(KNM_percent,KNM_avg)
KNM_comb <- KNM_comb[order(KNM_percent),]
KNM_comb <- data.frame(KNM_comb)
#Add KNM Points to Plot
KNM_plx <- predict(loess(KNM_comb$KNM_avg ~ KNM_comb$KNM_percent,
                         span = 0.3), se=T)
l <- msir::loess.sd(KNM_comb$KNM_avg ~ KNM_comb$KNM_percent, nsigma = 1.96)

#png("F07-duct-ccf-trends.png",width=12,height=6,units="in",res=1200)
#KM
plot(comb[,1], comb[,2], main = "Combined Duct Edge CCF Plots", 
     lwd = 2, ylim = c(0,1), xlab = "Edge Percent (%)", ylab = "CCF",pch=16)
lines(l$x, l$upper, lty=2,col='blue')
lines(l$x, l$lower, lty=2,col='blue')
lines(comb$percent,plx$fit,col='blue')
#KNM
points(KNM_comb[,1], KNM_comb[,2],col='gray',pch=16)
lines(l$x, l$upper, lty=2,col='red')
lines(l$x, l$lower, lty=2,col='red')
lines(KNM_comb$KNM_percent,KNM_plx$fit,col='red')

legend("bottomleft", title="Legend", c("KM","KM Trend","KM Prediction Band",
                                       "KNM","KNM Trend","KNM Prediction Band"), 
       col=c("black","blue","blue","gray","red","red"),
       lty = c(0,1,2,0,1,2),lwd = c(0,1,1,0,1,1),pch = c(16,NA,NA,16,NA,NA))
#dev.off()


#### ILLUSTRATE IMAGE PROCESSING/TRANSFORMATION ================================
#Figure 2 Components
setwd(duct_samples)
## STEP 1 - Read in Questioned Edge Image
samp1 <- load.image("KL.tiff")
# Rotate 90 degrees
samp1 <- imrotate(samp1, 90)
samp1 <- threshold(samp1, thr = "auto", approx = TRUE, adjust = 1.5)
# Edge Detection
edge1 <- cannyEdges(samp1)

# Plot with Identified Edge
plot(samp1)
highlight(edge1)

out1 <- data.frame(which(edge1, arr.ind = TRUE))
names(out1) <- c('x','y','a','b')
out1 <- out1 %>% 
  # Order data by x
  arrange(x) %>%
  # Subset the coord columns
  select(x, y)


## STEP 2 - Read in Known Edge Image
samp2 <- load.image("KR.tiff")
# Rotate 90 degrees
samp2 <- imrotate(samp2, 90)
samp2 <- threshold(samp2, thr = "auto", approx = TRUE, adjust = 1.5)
# Edge Detection
edge2 <- cannyEdges(samp2)

# Plot with Identified Edge - FIG 2(B)
plot(samp2)
highlight(edge2)

out2 <- data.frame(which(edge2, arr.ind = TRUE))
names(out2) <- c('x','y','a','b')
out2 <- out2 %>% 
  # Order data by x
  arrange(x) %>%
  # Subset the coord columns
  select(x, y)

## STEP 3 - Align at Xmin ----
if(out1$x[1] < out2$x[1]){
  shift <- out2$x[1] - out1$x[1]
  out1_align <- out1
  out2_align <- out2
  out2_align$x <- out2$x - shift
} else{
  shift <- out1$x[1] - out2$x[1]
  out1_align <- out1
  out2_align <- out2
  out1_align$x <- out1$x - shift
}

# Plot with Identified LM and Crop - FIG 2(C)
plot(out1_align, type = 'l')
lines(out2_align, col="red")

out1_crop <- out1_align[out1_align$x > min(out1_align$x) + 100 & 
                          out1_align$x < max(out1_align$x) - 100 ,]
out2_crop <- out2_align[out2_align$x > min(out2_align$x) + 100 & 
                          out2_align$x < max(out2_align$x) - 100 ,]
# plot(out1_crop, type = 'l',ylim = c(50, 300))
# lines(out2_crop, col='red')
out1_lm <- lm(y ~ x, data=out1_crop)
abline(lm(y ~ x, data=out1_crop))
out2_lm <- lm(y ~ x, data=out2_crop)
abline(lm(y ~ x, data=out2_crop), col="red")
abline(v = min(out1_align$x) + 100, col="blue") #crop line - left
abline(v = max(out1_align$x) - 100, col="blue") #crop line - right
legend("bottom", legend=c("KL", "KR", "Crop"),
       col=c("black", "red","blue"), lty=1, lwd = 2, cex=1)

out1_a <- out1_lm$coefficients[2]
out2_a <- out2_lm$coefficients[2]
out1_a - out2_a

out2_align <- data.frame(spdep::Rotation(out2, out1_a - out2_a))
names(out2_align) <- c('x','y')

out2_align$y <- out2_align$y + 200

#Figure 2(D)
#png("F02-alignmentD.png",width=12,height=6,units="in",res=1200)
plot(out1, type = 'l', main="TT1L and AAA1L", ylim = c(0,400))
lines(out2_align, col = "red")
legend("bottom", legend=c("AL1", "AR1"),
       col=c("black", "red"), lty=1, lwd = 2, cex=1)
dev.off()


#### CLASSIFIER ================================================================
dat <- do.call(rbind.data.frame, comb_out)
dat$type <- c()
pair <- data.frame(full=str_sub(dat$full, end=-2),
                   part=str_sub(dat$partial, end=-2))
for (i in 1:nrow(pair)) {
  if(pair$full[i]==pair$part[i]){
    dat$type[i] <- "KM"
  } else{dat$type[i] <- "KNM"}
}

dat$est_width <- floor(dat$partial_width * 100) + ceiling(dat$partial_width * 100) %% 2

KNM_widths <- filter(dat, type=="KNM")
KNM_widths <- unique(KNM_widths$est_width)
KM_widths <- filter(dat, type=="KM")
KM_widths <- unique(KM_widths$est_width)

result <- c()
pb = txtProgressBar(min = 0, max = nrow(dat), initial = 0, style = 3) 
for (i in 1:nrow(dat)) {
  setTxtProgressBar(pb,i)
  ccf_i <- dat$ccf[i]
  if(dat$type[i] == "KM"){
    KNMi <- filter(dat, type=="KNM",
                   est_width==dat$est_width[i])
    # In the event no cases of width exist, find nearest for comparison
    if(nrow(KNMi)==0){
      KNMi <- filter(dat, type=="KNM",
                    est_width==KNM_widths[which.min(abs(KNM_widths-dat$est_width[i]))])
    }
    if(ccf_i >= mean(KNMi$ccf)){
      result[i] <- "TP"
    } 
    if(ccf_i < mean(KNMi$ccf)){
      result[i] <- "FN"
    }
  } else{
    KMi <- filter(dat, type=="KM",
                  est_width==dat$est_width[i])
    # In the event no cases of width exist, find nearest for comparison
    if(nrow(KMi)==0){
      KMi <- filter(dat, type=="KM",
                    est_width==KM_widths[which.min(abs(KM_widths-dat$est_width[i]))])
    }
    if(ccf_i > mean(KMi$ccf)){
      result[i] <- "FP"
    } 
    if(ccf_i <= mean(KMi$ccf)){
      result[i] <- "TN"
    }
  }

}

head(dat)
load("dat_processed_full.RData")
#base::save.image("dat_processed_full.RData")
unique(dat$est_width)[order(unique(dat$est_width))]
dat_sub100 <- subset(dat, dat$est_width>=97)
dat_sub90 <- subset(dat, dat$est_width>=88 & dat$est_width<=92)
dat_sub70 <- subset(dat, dat$est_width>=68 & dat$est_width<=72)
dat_sub50 <- subset(dat, dat$est_width>=48 & dat$est_width<=52)
dat_sub30 <- subset(dat, dat$est_width>=28 & dat$est_width<=32)
dat_sub10 <- subset(dat, dat$est_width>=8 & dat$est_width<=12)

dat_sub <- dat_sub100
### Figure 4
png("km-knmdensity.png",width=9,height=8,units="in",res=1200)
ggplot(dat_sub100, aes(x=ccf, color=type,fill=type)) + 
  geom_density(alpha=.2) +
  theme(text = element_text(size=15)) +
  labs(x="CCF Score", y = "Density",fill="Type",color="Type")
dev.off()


### ERROR RATES ================================================================
TPR <- length(result[grep("TP",result)]) / (length(result[grep("TP",result)]) + 
                                              length(result[grep("FN",result)]))
FNR <- 1 - TPR
TNR <- length(result[grep("TN",result)]) / (length(result[grep("TN",result)]) + 
                                              length(result[grep("FP",result)]))
FPR <- 1 - TNR


library(ggplot2)
library("ggpubr")

#Figure 6
p10 <- ggplot(dat_sub10, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("10% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p10

p30 <- ggplot(dat_sub30, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("30% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p30

p50 <- ggplot(dat_sub50, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("50% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p50

p70 <- ggplot(dat_sub70, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("70% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p70

p90 <- ggplot(dat_sub90, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("90% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p90

p100 <- ggplot(dat_sub100, aes(x=ccf, fill=type)) +
  theme(text = element_text(size=10),legend.position = "none") +
  geom_density(alpha=0.4) +
  ggtitle("100% Edge") + labs(x = "CCF", y = "Density",fill="Type") + 
  scale_fill_grey() + theme_classic()
p100 #also Figure 4

### FIGURE 4
png("F04-density-plot.png",width=12,height=6,units="in",res=1200)
p100
dev.off()


### FIGURE 6
png("F06-edge-width-densities.png",width=12,height=6,units="in",res=1200)
ggarrange(p10,p30,p50,p70,p90,p100,
          ncol=3,nrow=2)
dev.off()


#### PLOT RMP USING ERROR RATES ================================================
library(scales)
library(grid)
library(gridExtra)

df <- data.frame(Width = unique(results$width)[1:20] * 100,
                 RMP = TPR)
df1 <- data.frame(RMP = round(TPR, digits = 3),
                  row.names = unique(results$width)[1:20] * 100)
df1 <- as.data.frame(t(df1))
plot2 <- ggplot(data = df, aes(x = Width, y = RMP)) + 
  geom_line(size = 1.5) +
  theme(text = element_text(size=20)) +
  # xlim(0,100) + 
  scale_x_continuous(name="Edge Width (%)", limits=c(0, 100), breaks = seq(0,100,by=10)) +
  labs(x="Edge Width (%)",y="Random Match Probability")
plot2

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1.1)),
  colhead = list(fg_params=list(cex = 1.1)),
  rowhead = list(fg_params=list(cex = 1.1)))

tbl <- tableGrob(df1, theme=mytheme)
png("RMP.png",width=13,height=7,units="in",res=300)
grid.arrange(plot2, tbl,
             nrow = 2, heights = c(2, 0.5))
dev.off()


#### DATA FILES ================================================================
#base::save.image("complete-data-plots.RData")
load("complete-data-plots.RData")
