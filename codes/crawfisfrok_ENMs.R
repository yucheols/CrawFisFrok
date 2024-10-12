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
options(java.parameters = "-Xmx8g")

# set random seed
set.seed(9)

# load packages we do {CLAP}
library(ENMeval)
library(megaSDM)
library(raster)
library(sf)
library(dplyr)
library(rasterVis)

#####  part 1 ::: get environmental data  ----------

# define clipping extent
ext <- c(-105.018,-84.675,28.311,44.003)

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


#####  part 3 ::: background data  ----------
bg <- dismo::randomPoints(mask = clim[[1]], n = 10000, p = cfis, excludep = T) %>% as.data.frame()
colnames(bg) = colnames(cfis)
head(bg)

# export bg
write.csv(bg, 'bg/bg.csv')


#####  part 4 ::: select environmental variables ----------

# make a correlation matrix
cor.mat <- raster::extract(clim, bg) %>% cor()
print(cor.mat) 

# run correlation test
ntbox::correlation_finder(cor_mat = cor.mat, threshold = 0.75, verbose = T)

# select low cor variables
clim.subs <- raster::stack(subset(clim, c('bio1', 'bio3', 'bio5', 'bio8', 'bio12', 'bio18')))
print(clim.subs)


#####  part 5 ::: test models  ----------
# "this pipeline is so charming" - our good friend Dr. JS

# test the candidate MaxEnt models 
cfis_mod <- ENMevaluate(taxon.name = 'crawfisfrok',
                        occs = cfis, 
                        envs = clim.subs, 
                        bg = bg, 
                        tune.args = list(fc = c('L','LQ','H','LQH','LQHP','LQHPT'), 
                                         rm = seq(1,5, by = 0.5)),
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

# find optimal model
cfis_opt <- cfis_res %>%
  dplyr::filter(or.10p.avg == min(or.10p.avg)) %>%
  dplyr::filter(auc.diff.avg == min(auc.diff.avg)) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg))

# check variable contribution
cfis_contrib <- eval.variable.importance(cfis_mod)[[cfis_opt$tune.args]]
print(cfis_contrib)

write.csv(cfis_contrib, 'outputs/contrib/crawfisfrok_optmodel_contrib.csv')

# look at the prediction map
cfis_opt_pred <- eval.predictions(cfis_mod)[[cfis_opt$tune.args]]
plot(cfis_opt_pred)

writeRaster(cfis_opt_pred, 'outputs/preds/crawfisfrok_optmodel_pred.tif')


#####  part 6 ::: response curves ----------

# resp curve data
cfis_resp_data <- ENMwrap::resp_data_pull(sp.name = 'crawfisfrok', model = eval.models(cfis_mod)[[cfis_opt$tune.args]], names.var = names(clim.subs))
head(cfis_resp_data)

# plot response
ENMwrap::plot_response(resp.data = cfis_resp_data)


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
cfis_lgm <- dismo::predict(object = eval.models(cfis_mod)[[cfis_opt$tune.args]], x = lgm, progress = 'text')
plot(cfis_lgm)


#####  part 8 ::: plot predictions ----------

# stack them
preds <- raster::stack(cfis_opt_pred, cfis_lgm)
names(preds) = c('Current', 'LGM')

# plot
gplot(preds) + 
  geom_tile(aes(fill = value)) + 
  coord_equal() + 
  facet_wrap(~ variable, ncol = 2, nrow = 1) + 
  scale_fill_gradientn(colors = rev(as.vector(pals::ocean.thermal(1000))), 
                       na.value = NA, 
                       name = "Suitability", 
                       breaks = c(0.1, 0.9), 
                       labels = c(paste0("Low: ", 0.1), paste0("High: ", 0.9))) + 
  xlab("Longitude") + ylab("Latitude") + 
  theme_bw() + 
  theme(strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 10)), 
        legend.text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.title.x = element_text(margin = margin(t = 15)), 
        axis.title.y = element_text(margin = margin(r = 15)), 
        axis.text = element_text(size = 12))

# export
ggsave('outputs/plots/preds_output.png', width = 20, height = 8, dpi = 600, units = 'cm')

