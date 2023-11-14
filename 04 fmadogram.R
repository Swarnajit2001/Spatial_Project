library(SpatialExtremes)
library(reshape2)
library(latex2exp)
library(patchwork)
load('We_Did_It.Rdata')


P = length(final.data)     #no. of indices
T = nrow(final.data[[1]])  #no. of years
N = ncol(final.data[[1]])  #no. of locations
#===========================================================================
#chi-measure for dependencies across the indices

chi1 = matrix(0, ncol = P, nrow = P)
time1 = matrix(0, ncol = P, nrow = P)

for(i in 1 : P){
  for(j in 1 : P){
    
    tock <- proc.time()[3]
    
    if(i == j){
      chi1[i,j] = 1
    }
    else if(i < j){
      
      
      
      
      H = numeric()
      
      for(l in 1 : ncol(final.data[[1]])){
  
        mat <- matrix(c(final.data[[i]][,l], final.data[[j]][,l]), ncol = 2)
        
        H[l] = unname((fmadogram(mat, coord = c(1,2)))[1,2])

      }
      
      chi.s = 2 - (1+2*H)/(1-2*H)

      chi1[i,j] = mean(chi.s)
    }
    else{
      
      chi1[i,j] = chi1[j,i]
      
    }
    print(paste(i,j))
    
    tick <- proc.time()[3]
    
    time1[i,j] = tick - tock
    
  }
}

name = c('Rx5day', 'R99p', 'Rx1day', 'R95p', 'R95pT', 'SDII', 'CWD', 'RPmm', 'PRCPTOT', 'R20mm')
colnames(chi1) = name
rownames(chi1) = name

save(chi1, time1, file = 'chi1.Rdata')

library(ggplot2)
library(reshape2)
final_chi1 = melt(chi1)

p1 = ggplot(final_chi1)+
  geom_tile(aes(x = Var1, y = Var2, fill = value))+
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 14))+
  scale_fill_gradient(high = '#021f38', low = '#a6d1f7')+
  labs(x ='',
       y = '',
       fill = '')

#===========================================================================
#chi-measure for dependencies across the space

# p = 1
chi2 = matrix(0, nrow = 27, ncol = P)

for(p in 1 : P){
  
  mat <- final.data[[p]]
  S = as.matrix(S.reg[,-3])
  
  mat.samp <- matrix(0, ncol = 1, nrow = nrow(mat))
  S.samp <- matrix(0, ncol = ncol(S), nrow = 1)
  for (r in 1:34){
    cur.ind <- which(S.reg[,3] == r)
    cur.ind.samp <- sample( cur.ind, ceiling(length(cur.ind) / 10) )
    mat.samp <- cbind(mat.samp, mat[ ,cur.ind.samp] ) 
    S.samp <- rbind(S.samp, S[cur.ind.samp, ] )
  }
  mat.samp <- mat.samp[,-1]
  S.samp <- S.samp[-1,]
  
  h = fmadogram(mat.samp,S.samp)
  h = na.omit(h)
  
  fmd = h[,2]
  foo = 2-(1+2*fmd)/(1-2*fmd)
  h[,2] = foo
  h = h[,-3]
  
  
  
  h.smooth <- matrix(0, nrow = 27, ncol = 2)
  for (i in 2:28)
  {
    ind <- which( ((h[,1] > i - 0.5) & (h[,1] < i + 0.5)) == TRUE )
    h.smooth[i-1,1] <- i
    h.smooth[i-1,2] <- mean( h[ind,2] )
  }
  
  chi2[,p] = h.smooth[,2]
  
  #----------------
  #aesthetics
  if (p <= 3){
    smile = ":-("
  }else if(p < 7){
    smile = ":-|"
  }else if(p < 10){
    smile = ":-)"
  }else{
    smile = ":-D"
  }
  print(paste0("We are at index ", p, ". ", smile))
}

name = c('Rx5day', 'R99p', 'Rx1day', 'R95p', 'R95pT', 'SDII', 'CWD', 'RPmm', 'PRCPTOT', 'R20mm')

dist = h.smooth[,1]
chi2 = cbind(dist, chi2)
colnames(chi2) = c('distance', name)
save(chi2, file = 'chi2.Rdata')

final_chi2 = melt(as.data.frame(chi2), id = 'distance')

library(ggplot2)
p2 = ggplot(final_chi2)+
  geom_line(aes(distance, value+0.1, col = variable), size = 1)+
  labs(x = 'h (degrees)',
       y = TeX(r"($\chi(h)$ )"),
       col = 'Indices')+
  xlim(0,30)+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17, hjust = 0.5),
        axis.title = element_text(size = 17))
p1+p2
  
