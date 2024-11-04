######  plot model output

# clear the working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(raster)
library(rgdal)
library(rasterVis)
library(ggplot2)


### load polygons
poly_cropped <- readOGR('shp/us_poly_cropped.shp')
us_states_crop <- readOGR('shp/us_states_cropped_tigris.shp')


###  load model predictions
cfis_oco_pred <- raster('outputs/preds/crawfisfrok_ocomodel_pred.tif')
cfis_lgm <- raster('outputs/preds/cfis_LGM_pred.tif')

preds <- stack(cfis_lgm, cfis_oco_pred)
names(preds) = c('LGM', 'Current')


###  plot
gplot(preds) + 
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable, ncol = 1, nrow = 2) + 
  scale_x_continuous(breaks = seq(-108, -79, by = 7)) +
  scale_fill_gradientn(colors = rev(as.vector(pals::ocean.thermal(1000))), 
                       na.value = NA, 
                       name = "Suitability", 
                       breaks = c(0.1, 0.9), 
                       labels = c(0.1, 0.9)) + 
  xlab("Longitude (°)") + ylab("Latitude (°)") + 
  geom_polygon(data = poly_cropped, aes(x = long, y = lat, group = group), color = 'darkgrey', linewidth = 0.5, fill = NA) +
  geom_polygon(data = us_states_crop, aes(x = long, y = lat, group = group), color = 'darkgrey', linewidth = 0.5, fill = NA) +
  theme_bw() + 
  theme(strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10, r = 10)), 
        legend.text = element_text(size = 12), 
        legend.position = 'top',
        axis.title = element_text(size = 14, face = "bold"), 
        axis.title.x = element_text(margin = margin(t = 15)), 
        axis.title.y = element_text(margin = margin(r = 15)), 
        axis.text = element_text(size = 12))

# export
ggsave('outputs/plots/preds_output_reorder.png', width = 10, height = 25, dpi = 600, units = 'cm')