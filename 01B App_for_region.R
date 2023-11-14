library(shiny)
require(rgdal)
require(ggplot2)
library(sp)

ui <- fluidPage(

    titlePanel("Choosing Regions of a Pixel"),

    sidebarLayout(
        sidebarPanel(
          
          numericInput("pixel", label = h3("Pixel input"), value = 1),
          hr()
        ),

        mainPanel(
           plotOutput("distPlot")
        )
    )
)


server <- function(input, output) {
  fn <- file.path("01A MeteorologicalSubDivisions.zip", fsep = "\\")
  utils::unzip(fn, exdir = tempdir())
  shp <- readOGR(dsn = file.path(tempdir(), ""), stringsAsFactors = F)
  
  #---------------------
  #map plot
  
  load("00 Data_Store_Raw Cleaned Data.Rdata")
  
  locs <- as.data.frame(S)
  #-------------
  #from lats longs to ne
  coordinates <- data.frame(lon = locs[,1], lat = locs[,2])
  spatial_points <- SpatialPoints(coordinates)
  proj4string(spatial_points) <- CRS("+proj=longlat +datum=WGS84")
  utm <- spTransform(spatial_points, CRS("+proj=utm +zone=43 +datum=WGS84"))
  northings <- coordinates(utm)[, 1]
  eastings <- coordinates(utm)[, 2]
  ret <- data.frame(nor = northings, est = eastings)
  locs.ne <- data.frame(nor = northings, est = eastings)
  #--------------
  locs.ne.1 = locs.ne
  locs.ne.1[,1] = locs.ne[,1]+3e6
  locs.ne.1[,2] = locs.ne[,2]+133e4

    output$distPlot <- renderPlot({
      
      cur.reg <- input$pixel
      
      ggplot() + geom_polygon(data = shp, aes(x = long, y = lat, group = group),
                              colour='black',fill=NA) + 
        geom_point(data = locs.ne.1, aes(x = nor, y = est), alpha = 0.3) + 
        geom_point(aes(x = locs.ne.1[cur.reg,1], y = locs.ne.1[cur.reg,2]), col = 'red') +
        ggtitle(cur.reg)
      
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
