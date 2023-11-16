library(ggmap)
library(ggplot2)
library(viridis)
library(Polychrome)

api_key = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
register_google(key = api_key)

bound <- c(
  left = min(S.reg$lon) - 4, bottom = min(S.reg$lat) - 4,
  right = max(S.reg$lon) + 4, top = max(S.reg$lat) + 4
)

india_centroid <- c(68.18625, 6, 97.41529, 37)

bdbox <- make_bbox(lon = c(67, 99), 
                   lat = c(39, 4))
india_map = get_map(location=bdbox, color="bw", zoom=5, maptype="terrain")


region_name = c('J&K','HP','PB','UK','HR','WR','WUP','ER','EUP',
                'BH','WB&S','A&M','ArP','NMMT','S&K','GJ','WMP','EMP','CG',
                'JK','GWB','K&G','MM','MT','VB','OR','TG','NIK','CK','SIK',
                'RS','CAP','KL','TN')
i = 1
S.reg[,3] = as.numeric(S.reg[,3])
for(i in 1 : nrow(S.reg)){
  S.reg[i,4] = region_name[S.reg[i,3]]
}


ggmap(india_map) + 
  geom_point(data = S.reg, aes(x = lon, y = lat, col = S.reg[,4]), size = 6, alpha = 0.6, shape = 15)+
  scale_fill_viridis(discrete = TRUE)+
  #scale_fill_manual(values = 1:34, name = region_name)+
  labs(x = 'Longitude',
       y = 'Latitude',
       col = 'Region')+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, hjust = 0.5))

