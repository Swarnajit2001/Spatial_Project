
#functions to get output of 10 variables

RX5day = function(data){
  return(max(sapply(1:361, function(i) mean(data[i:(i+4)])))) #applies the function on each element of vector
}

R99p = function(data, p99){
  return(sum(data[data>p99]))
}

Rx1day = function(data){return(max(data))}

R95p = function(data, p95){
  return(sum(data[data>p95]))
}

R95pT = function(data, p95){
  return(sum(data>p95))
}

SDII = function(data){
  return(sum(data)/sum(data>1))
} 

CWD = function(data){
  run_length = rle(data>1) #rle function counts the runlength of 1s and 0s
  return(max(run_length$lengths[run_length$values == 1]))
}

R10mm = function(data){
  return(sum(data>10))
}

PRCPTOT = function(data){
  return(sum(data[data>1]))
}

R20mm = function(data){
  return(sum(data>20))
}

#actual function
precip_vars = function(data, p99, p95){
  
  RX5day = function(data){
    return(max(sapply(1:361, function(i) mean(data[i:(i+4)])))) #applies the function on each element of vector
  }
  
  R99p = function(data, p99){
    return(sum(data[data>p99]))
  }
  
  Rx1day = function(data){return(max(data))}
  
  R95p = function(data, p95){
    return(sum(data[data>p95]))
  }
  
  R95pT = function(data, p95){
    return(sum(data>p95))
  }
  
  SDII = function(data){
    return(sum(data)/sum(data>1))
  } 
  
  CWD = function(data){
    run_length = rle(data>1) #rle function counts the runlength of 1s and 0s
    return(max(run_length$lengths[run_length$values == 1]))
  }
  
  R10mm = function(data){
    return(sum(data>10))
  }
  
  PRCPTOT = function(data){
    return(sum(data[data>1]))
  }
  
  R20mm = function(data){
    return(sum(data>20))
  }
  
  
  Output = list(RX5day = RX5day(data),
                R99p = R99p(data, p99),
                Rx1day = Rx1day(data),
                R95p = R95p(data, p95),
                R95pT = R95pT(data, p95),
                SDII = SDII(data),
                CWD = CWD(data),
                R10mm = R10mm(data),
                PRCPTOT = PRCPTOT(data),
                R20mm = R20mm(data))
  
  return(Output)
}

#===========================================================
#loading the data
library(tidyr)

data_0.25 = load('0.25grid_precip_data.Rdata')

all.data1 = list.Y

dat_vec = unlist(all.data1)
(p99 = quantile(dat_vec, .99, na.rm = T))
(p95 = quantile(dat_vec, .95, na.rm = T))

apply_on_mat = function(mat){
  a = apply(mat,1,FUN = function(x) precip_vars(x, p99, p95))
  b = unlist(a) %>% matrix(ncol = 10, byrow = T) #because we were getting a list inside every list, we are converting into a matrix as the output
  return(b)
}

# apply_on_mat(all.data1[[1]])

final_var_data = list(NULL)

for(i in 1 : length(all.data1)){
  final_var_data[[i]] = apply_on_mat(all.data1[[i]])
  Nan_index = which(final_var_data[[i]][,1] == 0)
  final_var_data[[i]][Nan_index,] = NA
  print(paste('Processed element', i))
}

save(final_var_data, file = 'precip_var_data.Rdata')


load('0.25grid_precip_data.Rdata')

library(readxl)
zone = as.character(read_excel('Regions.xlsx')$regions)

S.reg = data.frame(lon = S[,1], lat = S[,2], region = zone)


final.data <- list()

for (i in 1:ncol(final_var_data[[1]])){
  cur.mat <- matrix(0, nrow = length(final_var_data), ncol = nrow(final_var_data[[1]])) 
  
  for (j in 1:length(final_var_data)){
    cur.mat[j,] <- final_var_data[[j]][,i]
  }
  
  final.data[[i]] <- cur.mat
}

save(final.data, S.reg, file = 'We_Did_It.Rdata')
