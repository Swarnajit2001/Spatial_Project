library(EnvStats)
library(sn)
library(fGarch)
library(extraDistr)
library(evir)
library(ggplot2)

rm(list = ls())

load('03 We_Did_It.Rdata')

P = length(final.data)     #no. of indices
T = nrow(final.data[[1]])  #no. of years
N = ncol(final.data[[1]])  #no. of locations

#===========================================================================
#parameter estimations

final.par.est.nor <- matrix(0, ncol = 2, nrow = 10)
final.par.est.t <- matrix(0, ncol = 3, nrow = 10)
final.par.est.st <- matrix(0, ncol = 4, nrow = 10)
final.par.est.gev <- matrix(0, ncol = 3, nrow = 10)

for (i in 1:P)
{
  cur.par.est.nor <- matrix(0, ncol = 2, nrow = N)
  cur.par.est.t <- matrix(0, ncol = 3, nrow = N)
  cur.par.est.st <- matrix(0, ncol = 4, nrow = N)
  cur.par.est.gev <- matrix(0, ncol = 3, nrow = N)
  
  for (j in 1:N){
    
    cur.data <- final.data[[i]][,j]
    cur.data <- na.omit(cur.data)
    if (sum(which(cur.data == Inf)) > 0){
      cur.data <- cur.data[- which( cur.data == Inf )]
    }
    if (sum(which(cur.data == -Inf)) > 0){
      cur.data <- cur.data[- which( cur.data == -Inf )]
    }
    cur.data <- as.numeric(cur.data)
    
    
    #---
    #normal
    cur.par.est.nor[j, 1] <- mean(cur.data)
    n <- length(cur.data)
    cur.par.est.nor[j, 2] <- var(cur.data) * (n-1) / n
    
    #---
    #t
    cur.par.est.t[j, ] <-  unname(stdFit(cur.data)$par)
    
    #---
    #st
    cur.par.est.st[j, ] <- sstdFit(cur.data)$estimate
    
    #---
    #gev
    cur.par.est.gev[j, ] <- gev(cur.data)$par.ests
  
    #---
    #aesthetics
    print(paste(i,j)) 
  }
  
  cur.par.est.t <- cur.par.est.t[  cur.par.est.t[, 3] > 2, ]
  cur.par.est.st <- cur.par.est.st[  cur.par.est.st[, 3] > 2, ]
  cur.par.est.gev <- cur.par.est.gev[  abs(cur.par.est.gev[, 1]) < 0.5, ]
 
  
  
  final.par.est.nor[i, 1] <- mean( cur.par.est.nor[cur.par.est.nor[,1]<quantile(cur.par.est.nor[,1],0.2),1])
  final.par.est.nor[i, 2] <- mean( cur.par.est.nor[cur.par.est.nor[,2]<quantile(cur.par.est.nor[,2],0.2),2])
  
  final.par.est.t[i, 1] <- mean( cur.par.est.t[cur.par.est.t[,1]<quantile(cur.par.est.t[,1],0.3),1])
  final.par.est.t[i, 2] <- mean( cur.par.est.t[cur.par.est.t[,2]<quantile(cur.par.est.t[,2],0.3),2])
  final.par.est.t[i, 3] <- mean( cur.par.est.t[cur.par.est.t[,3]<quantile(cur.par.est.t[,3],0.3),3])
  
  final.par.est.st[i, 1] <- mean( cur.par.est.st[cur.par.est.st[,1]<quantile(cur.par.est.st[,1],0.3),1])
  final.par.est.st[i, 2] <- mean( cur.par.est.st[cur.par.est.st[,2]<quantile(cur.par.est.st[,2],0.3),2])
  final.par.est.st[i, 3] <- mean( cur.par.est.st[cur.par.est.st[,3]<quantile(cur.par.est.st[,3],0.3),3])
  final.par.est.st[i, 4] <- mean( cur.par.est.st[cur.par.est.st[,4]<quantile(cur.par.est.st[,4],0.3),4])
  
  final.par.est.gev[i, 1] <- mean( cur.par.est.gev[cur.par.est.gev[,1]<quantile(cur.par.est.gev[,1],0.3),1])
  final.par.est.gev[i, 2] <- mean( cur.par.est.gev[cur.par.est.gev[,2]<quantile(cur.par.est.gev[,2],0.3),2])
  final.par.est.gev[i, 3] <- mean( cur.par.est.gev[cur.par.est.gev[,3]<quantile(cur.par.est.gev[,3],0.3),3])
  
}

save(final.par.est.nor, final.par.est.t, final.par.est.st, final.par.est.gev, file = "Params_est.Rdata")

#==============================================================================
#plots

rm(list = ls())

load('03 We_Did_It.Rdata')
load("Params_est.Rdata")

P = length(final.data)

p = seq(0.001,0.999,0.001)

data.quantiles = function(mat){
  vec = as.vector(mat)
  vec =vec[vec<quantile(vec,0.3)]
  y = quantile(vec, p, na.rm = T)
  return(y)
}

y <- numeric()

for (i in 1:P){
  
  foo <- data.quantiles(final.data[[i]])
  foo <- unname(foo / max(foo))
  y <- c(y, foo)
}

#y = unname(unlist(lapply( final.data, FUN = data.quantiles) ) )

name = c('Rx5day', 'R99p', 'Rx1day', 'R95p', 'R95pT', 'SDII', 'CWD', 'RPmm', 'PRCPTOT','R20mm')
ind = rep(name, each = 999)

#---------
#Normal

quant1 = numeric()

for(i in 1 : P){
  foo <- qnorm(p, final.par.est.nor[i,1], sqrt(final.par.est.nor[i,2])) 
  foo <- foo / max(foo)
  quant1 = c(quant1, foo)
}

p1 <- ggplot()+
  geom_line(aes(quant1, y, col = ind), lwd = 0.6) + 
  geom_abline(slope=1, intercept=0, lwd = 0.8) + xlim(0,1) + 
  labs(x = "Theoritical quantiles", y = "Emprirical quantiles", col = "Indices") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 17),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom")

#---------
#t

quant2 = numeric()

for(i in 1 : P){
  foo <- qlst(p, final.par.est.t[i,3], final.par.est.t[i,1], final.par.est.t[i,2])
  foo <- foo / max(foo)
  quant2 = c(quant2, foo)
}

p2 <- ggplot()+
  geom_line(aes(quant2, y, col = ind), lwd = 0.6) + 
  geom_abline(slope=1, intercept=0, lwd = 0.8) + xlim(0,1) + 
  labs(x = "Theoritical quantiles", y = "Emprirical quantiles", col = "Indices") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 17),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom")

#---------
#Skew-t

library(sn)

quant3 = numeric()

for(i in 1 : P){
  foo <- qsstd(p, final.par.est.st[i,1], final.par.est.st[i,2],
               final.par.est.st[i,3], final.par.est.st[i,4])
  foo <- foo / max(foo)
  quant3 = c(quant3, foo)
}

p3 <- ggplot()+
  geom_line(aes(quant3, y, col = ind), lwd = 0.6) + 
  geom_abline(slope=1, intercept=0, lwd = 0.8) + xlim(0,1) + 
  labs(x = "Theoritical quantiles", y = "Emprirical quantiles", col = "Indices") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 17),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom")

#---------
#GEV

library(evd)

quant4 = numeric()

for(i in 1 : P){
  foo <- qgev(p, loc = final.par.est.gev[i,3], scale = final.par.est.gev[i,2], shape = final.par.est.gev[i,1])
  foo <- foo / max(foo)
  quant4 = c(quant4, foo)
}

p4 <- ggplot()+
  geom_line(aes(quant4, y, col = ind), lwd = 0.6) + 
  geom_abline(slope=1, intercept=0, lwd = 0.8) + xlim(0,1) + 
  labs(x = "Theoritical quantiles", y = "Emprirical quantiles", col = "Indices") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 17),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom")

#---------
#merging them together

library(patchwork)
library(gridExtra)
library(grid)

combined_plot.1 <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2) + 
  plot_annotation(theme = theme(legend.position = "top")) +
  plot_layout(guides = "collect") & xlab(NULL) & ylab(NULL)

combined_plot.1

combined_plot.2 <- wrap_elements(grid::textGrob("Emprirical quantiles",
                                                rot = 90, 
                                                x = 0.3, 
                                                y = 0.47, 
                                               gp = gpar(fontsize = 20))) + combined_plot.1  +
  plot_layout(ncol = 2, widths = c(0.1,2))

combined_plot.2

combined_plot.2 / wrap_elements(grid::textGrob("Theoritical quantiles", 
                                               x = 0.55, 
                                               y = 0.7, 
                                 gp = gpar(fontsize = 20))) +
  plot_layout(nrow = 2, heights = c(2,0.2))

i = 1 
 j = 1


#=================================================================
rm(list = ls()) 
 
load('03 We_Did_It.Rdata')
load("Params_est.Rdata")

P = length(final.data)     #no. of indices
T = nrow(final.data[[1]])  #no. of years
N = ncol(final.data[[1]])  #no. of locations
 
#only for skew-t
final.par.est.st.mean <- matrix(0, ncol = 4, nrow = 10)
final.par.est.st.cv <- matrix(0, ncol = 4, nrow = 10)

for (i in 1:P)
{
  cur.par.est.st <- matrix(0, ncol = 4, nrow = N)
  for (j in 1:N){
    cur.data <- final.data[[i]][,j]
    cur.data <- na.omit(cur.data)
    if (sum(which(cur.data == Inf)) > 0){
      cur.data <- cur.data[- which( cur.data == Inf )]
    }
    if (sum(which(cur.data == -Inf)) > 0){
      cur.data <- cur.data[- which( cur.data == -Inf )]
    }
    cur.data <- as.numeric(cur.data)
    #---
    #st
    cur.par.est.st[j, ] <- sstdFit(cur.data)$estimate
    cur.par.est.st[j,4] = cur.par.est.st[j,4]/sqrt(cur.par.est.st[j,2]) 
    
    #aesthetics
    print(paste(i,j)) 
    }
    cur.par.est.st <- cur.par.est.st[  cur.par.est.st[, 3] > 2, ]
    
    
    final.par.est.st.mean[i, 1] <- mean( cur.par.est.st[cur.par.est.st[,1]<quantile(cur.par.est.st[,1],0.3),1])
    final.par.est.st.mean[i, 2] <- mean( cur.par.est.st[cur.par.est.st[,2]<quantile(cur.par.est.st[,2],0.3),2])
    final.par.est.st.mean[i, 3] <- mean( cur.par.est.st[cur.par.est.st[,3]<quantile(cur.par.est.st[,3],0.3),3])
    final.par.est.st.mean[i, 4] <- mean( cur.par.est.st[cur.par.est.st[,4]<quantile(cur.par.est.st[,4],0.3),4])
    
    final.par.est.st.cv[i, 1] <- sd( cur.par.est.st[cur.par.est.st[,1]<quantile(cur.par.est.st[,1],0.3),1])
    final.par.est.st.cv[i, 2] <- sd( cur.par.est.st[cur.par.est.st[,2]<quantile(cur.par.est.st[,2],0.3),2])
    final.par.est.st.cv[i, 3] <- sd( cur.par.est.st[cur.par.est.st[,3]<quantile(cur.par.est.st[,3],0.3),3])
    final.par.est.st.cv[i, 4] <- sd( cur.par.est.st[cur.par.est.st[,4]<quantile(cur.par.est.st[,4],0.3),4])
    
    
}

final.par.est.st.cv = final.par.est.st.cv/final.par.est.st.mean 

save(final.par.est.st.mean, final.par.est.st.cv, file = 'only_for_skewt.Rdata')


final.par.est.st.mean
final.par.est.st.cv
