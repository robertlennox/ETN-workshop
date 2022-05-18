require(lubridate)
library(tidyverse)
library(osmdata) # package for working with streets
library(showtext) # for custom fonts
library(ggmap)
library(rvest)
library(sf)
library(data.table)
require(yaps)

# download data
# https://drive.google.com/drive/folders/1yiq68yp60Re8y4aIziXd3tsKsTG_1LWm?usp=sharing

tr <- opq(bbox = 'Aurland Norway') %>%
  add_osm_feature(key = 'name') %>%
  osmdata_sf()

aur<-tr$osm_polygons %>% 
  dplyr::filter(grepl("Vassb", name)) %>% 
  ggplot()+
  geom_sf(fill="#00aeff")+
  geom_sf(data=tr$osm_multipolygons %>% 
            dplyr::filter(grepl("Aurlandselv", name)),
          fill="#00aeff", colour="#00aeff")+
  geom_sf(data=tr$osm_multipolygons %>% 
            dplyr::filter(grepl("Aurlandsfjor", name)),
          fill="#00aeff")+
  theme_classic()+
  geom_sf(data=tr$osm_multipolygons %>% 
            dplyr::filter(grepl("selvi", name)),
          fill="#00aeff", colour=NA)+
  coord_sf(xlim=c(7.26, 7.31),
           ylim=c(60.862, 60.88))
setwd("C:/Users/robert.lennox/OneDrive - NINA/Aurland")
dets <- read.csv("august-yaps.csv") %>% # use the data I shared
  as_tibble %>% 
  dplyr::select(-1) %>% 
  mutate(dt=ymd_hms(dt))
hydros <- read.csv("C:/Users/robert.lennox/Downloads/hydros.csv") %>% 
  dplyr::select(-1)
temp <- read.csv("C:/Users/robert.lennox/OneDrive - NINA/Aurland/temp.csv") %>% 
  as_tibble %>% 
  mutate(dt=ymd_hms(dt))

# make the synchronization data frame: detections

detections<-dets %>% # sync tags
  dplyr::filter(ID==2055 |
                  ID==2057 |
                  ID==2058) %>% ## our sync tag was weird and had two IDs.. 
  mutate(ID=case_when(ID==as.numeric(2058)  ~ as.numeric(2057),
                      T~as.numeric(ID))) %>% 
  dplyr::select(dt, oid, Receiver, ID) %>% 
  mutate(epo=floor(as.numeric(dt))) %>% 
  dplyr::rename(ts=dt, tag=oid, serial=Receiver) %>% 
  dplyr::filter(month(ts)==8) %>%
  mutate(frac=seconds(ts)) %>% 
  mutate(frac2=as.character(frac)) %>% 
  separate(frac2, c("s", "fs")) %>% 
  mutate(fs=gsub("S", "", frac)) %>% 
  dplyr::select(-frac, -s) %>% 
  separate(fs, c("a", "b")) %>% 
  dplyr::select(-b) %>% 
  rename(frac=a) %>% 
  mutate(frac=paste0("0.", frac)) %>% 
  mutate(frac=as.numeric(frac)) %>% 
  mutate(frac=case_when(is.na(frac) ~ 0, T~frac)) %>%
  left_join(hydros %>% 
              mutate(serial2=serial) %>% 
              mutate(serial=as.integer(serial)) %>% 
              mutate(serial2=as.integer(serial2)) %>% 
              dplyr::select(serial, serial2)) %>% 
  dplyr::filter(!is.na(serial2)) %>% 
  dplyr::select(-tag) %>% 
  dplyr::rename(tag=ID)

detections %>% 
  left_join(hydros) %>% 
  group_by(x, tag) %>% 
  dplyr::summarise(n=n()) %>% 
  ggplot(aes(tag %>% factor, x %>% factor,  fill=n))+
  geom_tile()+
  scale_fill_viridis_c(option="plasma")+
  theme_classic()

## add in temperature data for speed of sound

r<-detections %>% 
  distinct(serial)
r<-r[[1]] # we just want a few receivers that are in the lake.. simple subset here

temp<-temp %>% 
  dplyr::filter(Receiver %in%  r) %>% 
  group_by(dt) %>% 
  dplyr::filter(year(dt)>2000) %>% 
  dplyr::summarise(temp=mean(temp)) 

temp<-tempToSs(temp=temp$temp, 
               sal=0, 
               depth=5) %>% 
  as_tibble %>% 
  bind_cols(temp) %>% 
  dplyr::rename(ss=value)

temp %>% 
  ggplot(aes(temp, ss))+
  theme_classic()+
  geom_line() # speed of sound and water temperature.. note that the high values are from the receiver being in air...

## make a syncrhonization list with hydrophones and sync tag detections

syncl<-list(hydros=hydros %>% setDT, # learn data.table
            detections=detections %>% 
              setDT)

max_epo_diff <- 200 # depends on burst interval of syncs... it is 600 secs
# trick is to align detection of each ping correctly
# should never be >50% of sync burst interval.. so never >300s
min_hydros <- 3 # toa matrix... if it doesn't work move to 2
time_keeper_idx <- 9
fixed_hydros_idx <- c(1:nrow(hydros))
n_offset_day <- 1 # with thelma or vemco, one usually enough.. maybe 2
n_ss_day <- 2 # ignored if not estimated
keep_rate<-1 # keep all the data.. might be ok with smaller dataset
# below one is a proportion of the data to keep
excl_self_detect = FALSE # I typically exclude own detections of sync tags, but since we only have two, we need them...
ss_data_what <- 'data'
ss_data<-temp %>% 
  dplyr::rename(ts=dt) %>% 
  setDT

inp_sync <- NULL
inp_sync <- getInpSync(syncl, 
                       max_epo_diff, 
                       min_hydros, 
                       time_keeper_idx, 
                       fixed_hydros_idx, 
                       n_offset_day, 
                       n_ss_day, 
                       ss_data=ss_data,
                       keep_rate=keep_rate, 
                       excl_self_detect=excl_self_detect, 
                       silent_check=F, 
                       ss_data_what=ss_data_what)

getSyncCoverage(inp_sync, plot=TRUE) # helps find timekeeper :)

sync_model_0 <- sync_model_1 <- sync_model_2 <- sync_model_3 <- sync_model_4 <- NULL
sync_model_0 <- getSyncModel(inp_sync, silent=TRUE, max_iter=1000, tmb_smartsearch = TRUE)
sync_model_1 <- fineTuneSyncModel(sync_model_0, eps_threshold=1E3, silent=TRUE)
sync_model_2 <- fineTuneSyncModel(sync_model_1, eps_threshold=1E2, silent=TRUE)
sync_model_3 <- fineTuneSyncModel(sync_model_2, eps_threshold=1E1, silent=TRUE)
sync_model_4 <- fineTuneSyncModel(sync_model_3, eps_threshold=1E1, silent=TRUE)
sync_model_5 <- fineTuneSyncModel(sync_model_4, eps_threshold=1E1, silent=TRUE)


setwd("C:/Users/robert.lennox/OneDrive - NINA/Aurland")
#sync_model<-readRDS("sync_Aurland2.rds")

sync_model<-sync_model_5

sync_model$eps_long %>% 
  as_tibble %>% 
  arrange(desc(E)) # check for outliers, rerun sync model with eps_threshold=1E1
# until the outliers are toast

sync_model$eps_long %>% 
  as_tibble %>%
  ggplot(aes(E, fill=hydro_idx %>% factor))+
  geom_histogram()+
  theme_classic()+
  labs(fill="Receiver")+
  theme(legend.position="top")

plotSyncModelResids(sync_model)+
  theme_classic()

sync_model$eps_long %>% 
  as_tibble %>% 
  ggplot(aes(factor(hydro_idx), E, fill=factor(sync_tag_idx)))+
  geom_boxplot()+
  theme_classic()+
  labs(fill="Sync Tag")+
  theme(legend.position="top")

aur+
  geom_point(data=sync_model$eps_long %>% 
               as_tibble %>% 
               dplyr::rename(idx=hydro_idx) %>% 
               left_join(hydros %>% 
                           st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
                           sf::st_transform(., crs = 4326) %>% 
                           as(., "Spatial") %>% 
                           as_tibble) %>% 
               group_by(coords.x1, coords.x2) %>% 
               dplyr::summarise(E=mean(E)),
             aes(coords.x1, coords.x2, size=E), colour="red")

##########
# Simple* YAPS
##########

fish_detections<-
  dets %>% 
  dplyr::filter(!is.na(oid)) %>% 
  dplyr::select(dt, oid, Receiver, ID) %>% 
  mutate(epo=floor(as.numeric(dt))) %>% 
  dplyr::rename(ts=dt, tag=oid, serial=Receiver) %>% 
  mutate(frac=seconds(ts)) %>% 
  mutate(frac2=as.character(frac)) %>% 
  separate(frac2, c("s", "fs")) %>% 
  mutate(fs=gsub("S", "", frac)) %>% 
  dplyr::select(-frac, -s) %>% 
  separate(fs, c("a", "b")) %>% 
  dplyr::select(-b) %>% 
  rename(frac=a) %>% 
  mutate(frac=paste0("0.", frac)) %>% 
  mutate(frac=as.numeric(frac)) %>% 
  mutate(frac=case_when(is.na(frac) ~ 0, T~frac)) %>% 
  dplyr::filter(!is.na(tag)) %>% 
  dplyr::filter(serial %in% hydros$serial)

dat <- applySync(toa=fish_detections %>%
                   dplyr::filter(yday(ts)==218) %>% 
                   dplyr::filter(tag==4712) %>% 
                   slice(1:500) %>% # only works with one tag at a time!
                   setDT, # never forget: Henrik <3 data.table
                 hydros=hydros, 
                 sync_model=sync_model)

sync_model$pl$OFFSET %>% 
  t() %>% 
  matplot()

dat %>% 
  as_tibble %>% 
  dplyr::select(ts, tag, serial, epo, eposync) %>% 
  ggplot(aes(eposync-epo))+
  geom_histogram(fill="black")+
  theme_classic()

matplot(toa[,7]-toa)


apply(toa,1,function(k) sum(!is.na(k))) %>% 
  plot()

hydros_yaps <- data.table::data.table(sync_model$pl$TRUE_H)
colnames(hydros_yaps) <- c('hx','hy','hz')

rbi_min <- 60
rbi_max <- 120

toa <- getToaYaps(synced_dat=dat, 
                  hydros=hydros_yaps, 
                  pingType='rbi', 
                  rbi_min=rbi_min, 
                  rbi_max=rbi_max)

nobs <- apply(toa, 1, function(k) sum(!is.na(k)))

aur+
  geom_point(data=toa %>% 
  as_tibble %>% 
  mutate(i=c(1:nrow(.))) %>% 
  gather(key, value, -i) %>% 
  dplyr::filter(!is.na(value)) %>% 
  dplyr::filter(i==5254) %>% 
  mutate(idx=str_sub(key, 2, -1)) %>% 
  mutate(idx=as.numeric(idx)) %>% 
  left_join(hydros %>% 
              st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
              sf::st_transform(., crs = 4326) %>% 
              as(., "Spatial") %>% 
              as_tibble),
aes(coords.x1, coords.x2, size=value-1633509498), colour="red")+
  labs(size="Time difference to detect")+
  theme(legend.position="top")

track<-runYaps(
  getInp(hydros_yaps, 
         toa, 
         E_dist="Mixture", # could be pure Gaussian or pure t.. dont mess
         n_ss=5,
         pingType="rbi", 
         sdInits=1, 
         rbi_min=rbi_min, 
         rbi_max=rbi_max, 
         ss_data_what="est", 
         ss_data=0, 
         bbox=100), # estimated positions cannot be 100m beyond
  silent=F, 
  tmb_smartsearch=T, 
  maxIter=5000)

# false convergence.. not ideal

# hey this does not look amazing, but we have gotten started at last

aur+
  geom_path(data=dat %>% 
               as_tibble %>% 
               left_join(hydros) %>% 
               st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
               sf::st_transform(., crs = 4326) %>% 
               as(., "Spatial") %>% 
               as_tibble %>% 
               group_by(ts=round_date(ts, "10 mins")) %>% 
               dplyr::summarise(coords.x1=mean(coords.x1),
                                coords.x2=mean(coords.x2)),
             aes(coords.x1, coords.x2), colour="grey")+
  geom_path(data=track$track %>% 
               as_tibble %>% 
               st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
               sf::st_transform(., crs = 4326) %>% 
               as(., "Spatial") %>% 
               as_tibble,
             aes(coords.x1, coords.x2))

## making YAPS fully operational
# the purrr::rerun function allows us to run the YAPS process five times
# we can extract the pseudo-AIC value and check the best fit, discard the rest

## that does not work because one bad run ruins the bunch.. enter:
# MAGIC YAPS

magicYAPS<-function(x){tryCatch({
  runYaps(
  getInp(hydros_yaps, 
         toa, 
         E_dist="Mixture", 
         n_ss=5, 
         pingType="rbi", 
         sdInits=1, 
         rbi_min=rbi_min, 
         rbi_max=rbi_max, 
         ss_data_what="est", 
         ss_data=0, 
         bbox=100), 
  silent=F, 
  tmb_smartsearch=TRUE, 
  maxIter=5000)},
  error=function(e){NA})}

YAPS_list<-5 %>% purrr::rerun(., magicYAPS())

YAPS_list %>% 
  purrr::map(purrr::pluck(8)) %>% 
  purrr::map(as_tibble) %>% 
  bind_rows(.id="id") 


YAPS_list %>% 
  purrr::map(purrr::pluck(8)) %>% 
  purrr::map(as_tibble) %>% 
  bind_rows(.id="id") %>% 
  ggplot(aes(top, x, colour=id))+
  geom_point()+
  theme_bw()

YAPS_list %>% 
  pluck(2) %>% 
  pluck(1) %>% 
  pluck(6)

aur+
  geom_path(data=YAPS_list %>% 
              purrr::map(purrr::pluck(8)) %>% 
              purrr::map(as_tibble) %>% 
              bind_rows(.id="id") %>% 
              st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
              sf::st_transform(., crs = 4326) %>% 
              as(., "Spatial") %>% 
              as_tibble,
            aes(coords.x1, coords.x2, colour=id),
            size=1.4)+
  scale_colour_manual(values=c("red", "orange", "yellow"))+
  theme(legend.position="top") 


## now let us see quantitatively what is the best track

# try my hacky code to get out the AIC* values for each YAPS run

fit<-YAPS_list %>%  # the magic number
  purrr::map(purrr::pluck(4)) %>% # get the AIC cols
  purrr::map(purrr::pluck(1)) %>% # take the number
  bind_cols() %>% # make a df
  t() %>% # oops wrong order
  as_tibble %>% 
  mutate(id=as.character(c(1, 3, 4))) %>% 
  dplyr::rename(AIC=V1)

aur+
  geom_path(data=YAPS_list %>% 
              purrr::map(purrr::pluck(8)) %>% 
              purrr::map(as_tibble) %>% 
              bind_rows(.id="id") %>% 
              st_as_sf(., coords = c("x", "y"), crs = 32633) %>% 
              sf::st_transform(., crs = 4326) %>% 
              as(., "Spatial") %>% 
              as_tibble %>% 
              left_join(fit),
            aes(coords.x1, coords.x2, colour=AIC %>% round %>% factor),
            size=1.4)+
  scale_colour_manual(values=c("red", "orange", "yellow"))+
  theme(legend.position="top") +
  labs(colour="AIC") # well.. that is discouraging, the best fit is obviously terrible



