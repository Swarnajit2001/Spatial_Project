#generated data

# data = rgamma(365,5,5)
# p99 = quantile(data, .99)
# p95 = quantile(data, .95)
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

precip_vars(data, p99, p95)
#================================================================================

#loading the data
library(tidyr)

main_data = load('Data_Store_Raw Cleaned Data.Rdata')

all.data1 = all.data

#removing columns of monthly and yearly totals (present in the data but not needed)

for(i in 1:length(all.data)){
  m = all.data[[i]]
  all.data1[[i]] = m[,-( (ncol(m)-1) : ncol(m) )]
}

#to determine the global 95th and 99th percentile\

dat_vec = unlist(all.data1)
(p99 = quantile(dat_vec, .99, na.rm = T))
(p95 = quantile(dat_vec, .95, na.rm = T))


#get the value of those 10 variables for all locations for all years
apply_on_mat = function(mat){
  a = apply(mat,1,FUN = function(x) precip_vars(x, p99, p95))
  b = unlist(a) %>% matrix(ncol = 10, byrow = T) #because we were getting a list inside every list, we are converting into a matrix as the output
  return(b)
}

final_var_data = list(NULL)

for(i in 1 : length(all.data1)){
  final_var_data[[i]] = apply_on_mat(all.data1[[i]])
  Nan_index = which(final_var_data[[i]][,1] == 0)
  final_var_data[[i]][Nan_index,] = NA
  print(paste('Processed element', i))
}

save(final_var_data, file = 'precip_var_data.Rdata')
# fna = 0
# for(i in 1 : 122){
#   fna[i] = sum(final_var_data[[i]] %>% rowSums() == 0)
# }
# fna


