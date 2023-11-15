rm(list = ls())

library(fda)

load("03 We_Did_It.Rdata")

#----------------------------------------------------------------------------

N.all <- ncol(final.data[[1]])

P <- 10
all.reg <- c(1:5, 7:8, 10:19, 21:28, 30, 31, 33, 34)

final.mat <- matrix(0, ncol = P, nrow = N.all)

#creating x.t's
L <- 5
basis <-  create.bspline.basis(c(1901,2022), nbasis = L)
B.mat <- eval.basis(1901:2022, basis)
X <- t(B.mat)

x.1901 <- X[,1]
x.2022 <- X[,ncol(X)]

for (i in all.reg)
{

  S.indv = which(S.reg[,3] == i) 
  
  N = length(S.indv)
  
  reg.name <- ifelse(i < 9, paste0("0",i), as.character(i))
  file.name <- paste0("final.MCMC.R", reg.name, ".Rdata")
  load(file.name)

  beta.chain <- out$beta[2501:5000, ]
  beta.mean <- colMeans(beta.chain)
  #beta.var <- apply(beta.chain, 2, var)
  
  X.1901 <- diag(N * P) %x% x.1901
  X.2022 <- diag(N * P) %x% x.2022
  
  mu.2022 <- t(X.2022) %*% beta.mean
  mu.1901 <- t(X.1901) %*% beta.mean
  
  mu.mat.2022 <- matrix(mu.2022, ncol = N)
  mu.mat.1901 <- matrix(mu.1901, ncol = N)
  
  change.mat <- (mu.mat.2022 - mu.mat.1901) / 500
  
  final.mat[S.indv, ] <- t(change.mat)
  
  rm(out)
  
  #aesthetics
  print(paste0("We are at region = ", i))

}

save(final.mat, file = "07 MCMC_Plots.Rdata")

#=====================================================================


# register_google(key = api_key)
# 
# bound <- c(
#   left = min(S.reg$lon) - 4, bottom = min(S.reg$lat) - 4,
#   right = max(S.reg$lon) + 4, top = max(S.reg$lat) + 4
# )
# 
# india_centroid <- c(68.18625, 6, 97.41529, 37)
# 
# bdbox <- make_bbox(lon = c(67, 99), 
#                    lat = c(39, 4))
# india_map = get_map(location=bdbox, color="bw", zoom=5, maptype="terrain")
# 
# save(india_map, file = "India_map.Rdata")

#-- -- -- -- -- -- -- -- -- --

rm(list = ls())

library(ggmap)
library(ggplot2)
library(viridis)
library(Polychrome)

load("07 India_map.Rdata")
load("07 MCMC_Plots.Rdata")
load("03 We_Did_It.Rdata")

coords <- S.reg[,1:2]

#-------------
#var.1:

#df.01 <- cbind(coords, scale(final.mat[,1], center = T,scale = T))
df.01 <- cbind(coords, final.mat[,1])

colnames(df.01) <- c("lon", "lat", "var")

df.01[which(df.01$var > quantile(df.01$var, 0.75)), 3] = quantile(df.01$var, 0.75)
df.01[which(df.01$var < quantile(df.01$var, 0.25)), 3] = quantile(df.01$var, 0.25)

p1 <- ggmap(india_map) + 
  geom_point(data = df.01, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'Rx5day') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p1

#-------------
#var.2:

df.02 <- cbind(coords, final.mat[,2])

colnames(df.02) <- c("lon", "lat", "var")

df.02[which(df.02$var > quantile(df.02$var, 0.75)), 3] = quantile(df.02$var, 0.75)
df.02[which(df.02$var < quantile(df.02$var, 0.25)), 3] = quantile(df.02$var, 0.25)

p2 <- ggmap(india_map) + 
  geom_point(data = df.02, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'R99p') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))


#p2

#-------------
#var.3:

df.03 <- cbind(coords, final.mat[,3])

colnames(df.03) <- c("lon", "lat", "var")

df.03[which(df.03$var > quantile(df.03$var, 0.75)), 3] = quantile(df.03$var, 0.75)
df.03[which(df.03$var < quantile(df.03$var, 0.25)), 3] = quantile(df.03$var, 0.25)

p3 <- ggmap(india_map) + 
  geom_point(data = df.03, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'Rx1day') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))
#p3

#-------------
#var.4:

df.04 <- cbind(coords, final.mat[,4])

colnames(df.04) <- c("lon", "lat", "var")

df.04[which(df.04$var > quantile(df.04$var, 0.75)), 3] = quantile(df.04$var, 0.75)
df.04[which(df.04$var < quantile(df.04$var, 0.25)), 3] = quantile(df.04$var, 0.25)

p4 <- ggmap(india_map) + 
  geom_point(data = df.04, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'R95p') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p4

#-------------
#var.5:

df.05 <- cbind(coords, final.mat[,5])

colnames(df.05) <- c("lon", "lat", "var")

df.05[which(df.05$var > quantile(df.05$var, 0.75)), 3] = quantile(df.05$var, 0.75)
df.05[which(df.05$var < quantile(df.05$var, 0.25)), 3] = quantile(df.05$var, 0.25)

p5 <- ggmap(india_map) + 
  geom_point(data = df.05, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'R95pT') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p5

#-------------
#var.6:

df.06 <- cbind(coords, final.mat[,6])

colnames(df.06) <- c("lon", "lat", "var")

df.06[which(df.06$var > quantile(df.06$var, 0.75)), 3] = quantile(df.06$var, 0.75)
df.06[which(df.06$var < quantile(df.06$var, 0.25)), 3] = quantile(df.06$var, 0.25)

p6 <- ggmap(india_map) + 
  geom_point(data = df.06, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'SDII') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p6

#-------------
#var.7:

df.07 <- cbind(coords, final.mat[,7])

colnames(df.07) <- c("lon", "lat", "var")

df.07[which(df.07$var > quantile(df.07$var, 0.75)), 3] = quantile(df.07$var, 0.75)
df.07[which(df.07$var < quantile(df.07$var, 0.25)), 3] = quantile(df.07$var, 0.25)

p7 <- ggmap(india_map) + 
  geom_point(data = df.07, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'CWD') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p7

#-------------
#var.8:

df.08 <- cbind(coords, final.mat[,8])

colnames(df.08) <- c("lon", "lat", "var")

df.08[which(df.08$var > quantile(df.08$var, 0.75)), 3] = quantile(df.08$var, 0.75)
df.08[which(df.08$var < quantile(df.08$var, 0.25)), 3] = quantile(df.08$var, 0.25)

p8 <- ggmap(india_map) + 
  geom_point(data = df.08, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'RPmm') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p8

#-------------
#var.9:

df.09 <- cbind(coords, final.mat[,9])

colnames(df.09) <- c("lon", "lat", "var")

df.09[which(df.09$var > quantile(df.09$var, 0.75)), 3] = quantile(df.09$var, 0.75)
df.09[which(df.09$var < quantile(df.09$var, 0.25)), 3] = quantile(df.09$var, 0.25)

p9 <- ggmap(india_map) + 
  geom_point(data = df.09, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'PRCPTOT') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p9

#-------------
#var.10:

df.10 <- cbind(coords, final.mat[,10])

colnames(df.10) <- c("lon", "lat", "var")

df.10[which(df.10$var > quantile(df.10$var, 0.75)), 3] = quantile(df.10$var, 0.75)
df.10[which(df.10$var < quantile(df.10$var, 0.25)), 3] = quantile(df.10$var, 0.25)

p10 <- ggmap(india_map) + 
  geom_point(data = df.10, mapping = aes(x = lon, y = lat, col = var), size = 2, alpha = 0.6, shape = 15) +
  scale_fill_viridis(discrete = TRUE)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'R20mm') +
  scale_color_viridis() +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, hjust = 0.5))

#p10

#==========================================================================
#all together

library(patchwork)
library(gridExtra)
library(grid)

combined_plot.1 <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) + 
  plot_annotation(theme = theme(legend.position = "top"))

combined_plot.1

combined_plot.2 <- p7 + p8 + p9 + p10 + plot_layout(ncol = 2) + 
  plot_annotation(theme = theme(legend.position = "top"))

combined_plot.2

