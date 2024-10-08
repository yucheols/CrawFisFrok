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

# load packages we do {CLAP}
library(ENMwrap)
library(raster)
library(sf)
library(dplyr)
library(megaSDM)


#####  part 1 ::: get environmental data  ----------

# define clipping extent
ext <- c(-105.018,-84.675,28.311,44.003)

# WorldClim 5km raster
clim <- raster::stack(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
clim <- raster::crop(clim, ext)

plot(clim[[1]])


#####  part 2 ::: get occurrence points  ----------

# collect occurrence points 
OccurrenceCollection(spplist = c('Lithobates areolatus',
                                 'Lithobates areolatus areolatus',
                                 'Lithobates areolatus circulosus'),
                     output = 'occs',
                     trainingarea = extent(clim[[1]]))
