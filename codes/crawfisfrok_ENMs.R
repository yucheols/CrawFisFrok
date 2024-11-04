##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##                                              __  _.---""---._  __                                        ##                              
##                                          /####\/              \/####\              #######               ##
##                                         (   /\ )              ( /\   )             ##                    ##
##                                          \____/                \____/              #####                 ##
##                                        _/                          \_              ##                    ##
##                                    .-"    .                      .    "-.          ##                    ##
##                                    |  |   \.._                _../   |  |                                ##
##                                     \  \    \."-.__________.-"./    /  /                                 ##
##                                        \  \    "--.________.--"    /  /            #####                 ##
##                                      __\  \                    /  /__              #   ##                ##
##                                   ./    )))))                  (((((    \.         #####                 ##  
##                                    \                                      /        #   #                 ##
##                                      \           \_          _/           /        #    #                ##
##                                       \    \____/""-.____.-""\____/    /                                 ##
##                                         \    \                  /    /                                   ##
##                                         \.  .|                |.  ./                ####                 ##
##                                        ." / |  \              /  | \  ".           #    #                ##
##                                      ."  /   |   \            /   |   \   ".       #    #                ##
##                                      /.-./.--.|.--.\          /.--.|.--.\.-.|       ####                 ##
##                                                                                                          ##
##                                                                                                          ##
##                                                                                     ##   ##              ## 
##                                                                                     ## ##                ##
##                                                                                     ##                   ##   
##                                                                                     ###                  ##
##                                                                                     ## ##                ## 
##                                                                                     ##   ##              ##
##                                                                                                          ##
##                                                                                                          ##
##############################################################################################################
##############################################################################################################


###########    CrawFisFrok ENM == We go NOW

# kindly clear the working environment
rm(list = ls(all.names = T))
gc()

# use disk-based processing to avoid loading everything into memory at once
rasterOptions(tmpdir = 'tmpdir', progress = 'text', chunksize = 1e+07, maxmemory = 1e+09)

# increase Java heap space
options(java.parameters = "-Xmx10g")

# set random seed
set.seed(9)

# load packages we do {CLAP}
library(ENMeval)
library(ENMwrap)
library(megaSDM)
library(raster)
library(sf)
library(dplyr)
library(rasterVis)
library(ggplot2)
library(rnaturalearth)
library(tigris)

#####  part 1 ::: get environmental data  ----------

# define clipping extent
ext <- c(-108.8371,-78.6800,19.0456,50.2829)

# WorldClim 5km raster
clim <- raster::stack(list.files(path = 'E:/env layers/CHELSA_cur_V1_2B_r30s/30sec/', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, ext)

names(clim) = gsub('_', '', names(clim))
print(clim)

plot(clim[[1]])

# export processed layers
for (i in 1:nlayers(clim)) {
  r <- clim[[i]]
  raster::writeRaster(r, paste0('clim_processed/', names(clim)[i], '.bil'), overwrite = T)
}

# import shortcut
clim <- raster::stack(list.files(path = 'clim_processed/', pattern = '.bil$', full.names = T))


#####  part 2 ::: get occurrence points  ----------

# collect occurrence points == downloaded on 9 Oct 2024 at 2:25 PM
#OccurrenceCollection(spplist = c('Lithobates areolatus',
#                                 'Lithobates areolatus areolatus',
#                                 'Lithobates areolatus circulosus'),
#                     output = 'occs',
#                     trainingarea = extent(clim[[1]]))

# load occurrence points
cfis <- read.csv('occs/Lithobates_areolatus.csv') %>% select('species', 'decimalLongitude', 'decimalLatitude')
colnames(cfis) = c('species', 'long', 'lat')
head(cfis)

# thin occurrence points
cfis <- SDMtune::thinData(coords = cfis[, c(2,3)], env = terra::rast(clim[[1]]), x = 'long', y = 'lat', verbose = T, progress = T)
head(cfis)

points(cfis)
nrow(cfis)

# export thinned
write.csv(cfis, 'occs/cfis_thinned.csv')


#####  part 3 ::: background data  ----------

# generate buffer == 300 km around the occurrence points
buff <- ENMwrap::buff_maker(occs_list = list(cfis), envs = clim[[1]], buff_dist = 300000)
plot(buff[[1]], border = 'blue', lwd = 2, add = T)

# sample bg
bg <- ENMwrap::bg_sampler(envs = clim[[1]], n = 10000, occs_list = list(cfis), buffer_list = buff, method = 'buffer') %>% dplyr::bind_rows()
colnames(bg) = colnames(cfis)
head(bg)

# export bg
write.csv(bg, 'bg/bg.csv')


#####  part 4 ::: select environmental variables ----------

# make a correlation matrix
cor.mat <- raster::extract(clim, bg) %>% cor()
print(cor.mat) 

# run correlation test
ntbox::correlation_finder(cor_mat = cor.mat, threshold = 0.70, verbose = T)

# select low cor variables
clim.subs <- raster::stack(subset(clim, c('bio1','bio3','bio5','bio8','bio12','bio18')))
print(clim.subs)


#####  part 5 ::: test models  ----------
# "this pipeline is so charming" - our good friend Dr. JS

# test the candidate MaxEnt models 
cfis_mod <- ENMevaluate(taxon.name = 'crawfisfrok',
                        occs = cfis, 
                        envs = clim.subs, 
                        bg = bg, 
                        tune.args = list(fc = c('L','LQ','H','LQH','LQHP','LQHPT'), 
                                         rm = seq(1, 4, by = 0.5)),
                        partitions = 'block',
                        partition.settings = list(orientation = 'lat_lon'),
                        algorithm = 'maxent.jar',
                        doClamp = T,
                        parallel = T,
                        parallelType = 'doSNOW')

# save model object
saveRDS(cfis_mod, 'outputs/models/crawfisfrok_tuning.rds')

# check the model results
cfis_res <- eval.results(cfis_mod)
print(cfis_res)

write.csv(cfis_res, 'outputs/models/crawfisfrok_tuning_results.csv')

# find optimal model == which model is o'co??
cfis_oco <- cfis_res %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(auc.diff.avg == min(auc.diff.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg)) %>%
  print()

# check variable contribution
cfis_contrib <- eval.variable.importance(cfis_mod)[[cfis_oco$tune.args]]
print(cfis_contrib)

write.csv(cfis_contrib, 'outputs/contrib/crawfisfrok_ocomodel_contrib.csv')

# look at the prediction map
cfis_oco_pred <- eval.predictions(cfis_mod)[[cfis_oco$tune.args]]
plot(cfis_oco_pred)

writeRaster(cfis_oco_pred, 'outputs/preds/crawfisfrok_ocomodel_pred.tif', overwrite = T)


#####  part 6 ::: response curves ----------

# resp curve data
cfis_resp_data <- ENMwrap::resp_data_pull(sp.name = 'crawfisfrok', model = eval.models(cfis_mod)[[cfis_oco$tune.args]], names.var = names(clim.subs))
head(cfis_resp_data)

# reorder plotting order
cfis_resp_data$var <- factor(cfis_resp_data$var, levels = c('bio1', 'bio3', 'bio5', 'bio8', 'bio12', 'bio18'))

# plot response
ENMwrap::plot_response(resp.data = cfis_resp_data) +
  theme(axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 16, face = 'italic'),
        legend.position = 'top')

# export response curves
ggsave('outputs/plots/cfis_resp_curves.png', width = 30, height = 22, dpi = 800, units = 'cm')


#####  part 7 ::: hindcasting ----------

# import LGM layers
lgm <- raster::stack(list.files(path = 'E:/env layers/chelsa_LGM_v1_2B_r30s/30sec', pattern = '.tif$', full.names = T))
lgm <- raster::crop(lgm, ext)
plot(lgm[[1]])

names(lgm) = gsub('_', '', names(lgm))

# subset
lgm <- raster::stack(subset(lgm, names(clim.subs)))
print(lgm)

# check climate units
ENMwrap::unit_check(ref.env = clim.subs, proj.env = lgm, n = 10000)

# it o'co proceed to projection
cfis_lgm <- dismo::predict(object = eval.models(cfis_mod)[[cfis_oco$tune.args]], x = lgm, progress = 'text')
plot(cfis_lgm)

# export LGM prediction
writeRaster(cfis_lgm, 'outputs/preds/cfis_LGM_pred.tif', overwrite = T)


#####  part 8 ::: plot predictions ----------

###  prep global polygon data
# get polygon data
glob_poly <- ne_countries(scale = 'small', returnclass = 'sf')
glob_poly <- st_transform(glob_poly, crs = crs(clim.subs))

# check geometry
st_is_valid(glob_poly, reason = T)

# fix invalid geometry
sf_use_s2(F)
glob_poly <- st_make_valid(glob_poly)

# crop to the modeling extent
poly_bbox <- st_bbox(extent(clim.subs), crs = st_crs(glob_poly))
poly_cropped <- st_crop(glob_poly, poly_bbox)


### prep US polygon data
# get data
us_states <- states()
us_states <- st_transform(us_states, crs = crs(clim.subs))

# crop
us_states_crop <- st_crop(us_states, poly_bbox)


###  plot
# stack preds
preds <- raster::stack(cfis_lgm, cfis_oco_pred)
names(preds) = c('LGM', 'Current')

# plot
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
  xlab("Longitude") + ylab("Latitude") + 
  geom_sf(data = poly_cropped, inherit.aes = F, fill = NA, color = 'darkgrey', lwd = 0.5) +
  geom_sf(data = us_states_crop, inherit.aes = F, fill = NA, color = 'darkgrey', lwd = 0.5) +
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
ggsave('outputs/plots/preds_output.png', width = 10, height = 25, dpi = 600, units = 'cm')


#####  part 8 ::: MESS ----------

# run MESS
cfis_mess <- ntbox::ntb_mess(M_stack = clim.subs, G_stack = lgm)
print(cfis_mess)

# plot MESS
gplot(cfis_mess) +
  geom_tile(aes(fill = value)) + 
  coord_equal() + 
  scale_fill_gradientn(colors = as.vector(pals::coolwarm(1000)), 
                       na.value = NA,
                       breaks = c(-75, 60),
                       labels = c('High', 'Low'),
                       name = "MESS",
                       trans = 'reverse') +
  xlab("Longitude") + ylab("Latitude") +
  geom_sf(data = poly_cropped, inherit.aes = F, fill = NA, color = 'darkgrey', lwd = 0.5) +
  geom_sf(data = us_states_crop, inherit.aes = F, fill = NA, color = 'darkgrey', lwd = 0.5) +
  theme_bw()+
  theme(strip.text = element_text(size = 14, face = "italic"), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.title.x = element_text(margin = margin(t = 15)), 
        axis.title.y = element_text(margin = margin(r = 15)), 
        axis.text = element_text(size = 12))

# export the plot
ggsave('outputs/plots/MESS.png', width = 20, height = 18, dpi = 600, units = 'cm')

# export the MESS layer
writeRaster(cfis_mess, 'outputs/MESS/cfis_MESS.tif', overwrite = T)
