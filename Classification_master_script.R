# baseline through hanna meyer method
library(tidymodels)
library(sf)
library(mapview)
library(raster)
library(here)
library(fasterize)
library(randomForest)
library(doParallel)
library(parallel)
library(CAST)
library(caret)
library(exactextractr)
library(plotly)
library(stringr)

#### Setup band names, location, and period ####

names_S2 <- c('red', 'green', 'blue', 'nir', 'swir_1', 'swir_2', 'red_edge_1', 'red_edge_2', "smad", "emad", "bcmad", "NDVI", "BSI", "MNDWI", "CIRE")
names_S1 <- c("VV", "VH", 'RVI')

location <- "Catandica"

#### Load rasters ####

##### 1x12 months #####
mo_12 <- brick(here("data","DEA_round2", paste0(location,"_mo12_20192020.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_12 <- (mo_12[[4]]/mo_12[[8]]) - 1
mo_12 <- addLayer(mo_12, cire_12 ) 
names(mo_12) <- paste0(names_S2)


mo_12_SAR <- brick(here("data","DEA_round2",paste0(location,"_12mo_20192020_SAR.tif")))
names(mo_12_SAR) <- paste0(names_S1)

##### 2x6 months #####
mo_6_Q41 <- brick(here("data","DEA_round2", paste0(location,"_6mo_2019Q4_2020Q1.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_6_Q41 <- (mo_6_Q41[[4]]/mo_6_Q41[[8]]) - 1
mo_6_Q41 <- addLayer(mo_6_Q41, cire_6_Q41 ) 
names(mo_6_Q41) <- paste0(names_S2, "_6m_Q41")

mo_6_Q23 <- brick(here("data","DEA_round2", paste0(location,"_6mo_2020Q2_2020Q3.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_6_Q23 <- (mo_6_Q23[[4]]/mo_6_Q23[[8]]) - 1
mo_6_Q23 <- addLayer(mo_6_Q23, cire_6_Q23 ) 
names(mo_6_Q23) <- paste0(names_S2, "_6m_Q23")


mo_6_Q41_SAR <- brick(here("data","DEA_round2",paste0(location,"_6mo_2019Q4_2020Q1_SAR.tif")))
names(mo_6_Q41_SAR) <- paste0(names_S1, "_6m_Q41")

mo_6_Q23_SAR <- brick(here("data","DEA_round2",paste0(location,"_6mo_2020Q2_2020Q3_SAR.tif")))
names(mo_6_Q23_SAR) <- paste0(names_S1, "_6m_Q23")

##### 4x3 months #####
mo_3_Q4 <- brick(here("data","DEA_round2", paste0(location,"_3mo_2019_Q4.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_3_Q4 <- (mo_3_Q4[[4]]/mo_3_Q4[[8]]) - 1
mo_3_Q4 <- addLayer(mo_3_Q4, cire_3_Q4 ) 
names(mo_3_Q4) <- paste0(names_S2, "_3m_Q4")

mo_3_Q1 <- brick(here("data","DEA_round2", paste0(location,"_3mo_2020_Q1.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_3_Q1 <- (mo_3_Q1[[4]]/mo_3_Q1[[8]]) - 1
mo_3_Q1 <- addLayer(mo_3_Q1, cire_3_Q1) 
names(mo_3_Q1) <- paste0(names_S2, "_3m_Q1")

mo_3_Q2 <- brick(here("data","DEA_round2", paste0(location,"_3mo_2020_Q2.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_3_Q2 <- (mo_3_Q2[[4]]/mo_3_Q2[[8]]) - 1
mo_3_Q2 <- addLayer(mo_3_Q2, cire_3_Q2 ) 
names(mo_3_Q2) <- paste0(names_S2, "_3m_Q2")

mo_3_Q3 <- brick(here("data","DEA_round2", paste0(location,"_3mo_2020_Q3.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_3_Q3 <- (mo_3_Q3[[4]]/mo_3_Q3[[8]]) - 1
mo_3_Q3 <- addLayer(mo_3_Q3, cire_3_Q3 ) 
names(mo_3_Q3) <- paste0(names_S2, "_3m_Q3")


mo_3_Q4_SAR <- brick(here("data","DEA_round2",paste0(location,"_3mo_2019_Q4_SAR.tif")))
names(mo_3_Q4_SAR) <- paste0(names_S1, "_3m_Q4")

mo_3_Q1_SAR <- brick(here("data","DEA_round2",paste0(location,"_3mo_2020_Q1_SAR.tif")))
names(mo_3_Q1_SAR) <- paste0(names_S1, "_3m_Q1")

mo_3_Q2_SAR <- brick(here("data","DEA_round2",paste0(location,"_3mo_2020_Q2_SAR.tif")))
names(mo_3_Q2_SAR) <- paste0(names_S1, "_3m_Q2")

mo_3_Q3_SAR <- brick(here("data","DEA_round2",paste0(location,"_3mo_2020_Q3_SAR.tif")))
names(mo_3_Q3_SAR) <- paste0(names_S1, "_3m_Q3")

##### 6x2 months #####
mo_2_2019_10_11 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2019_10_11.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2019_10_11 <- (mo_2_2019_10_11[[4]]/mo_2_2019_10_11[[8]]) - 1
mo_2_2019_10_11 <- addLayer(mo_2_2019_10_11, cire_2_2019_10_11 ) 
names(mo_2_2019_10_11) <- paste0(names_S2, "_2m_2019_10_11")

mo_2_2020_12_01 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2020_12_01.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2020_01_02 <- (mo_2_2020_12_01[[4]]/mo_2_2020_12_01[[8]]) - 1
mo_2_2020_12_01 <- addLayer(mo_2_2020_12_01, cire_2_2020_01_02 ) 
names(mo_2_2020_12_01) <- paste0(names_S2, "_2mo_2020_12_01")

mo_2_2020_02_03 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2020_02_03.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2020_03_04 <- (mo_2_2020_02_03[[4]]/mo_2_2020_02_03[[8]]) - 1
mo_2_2020_02_03 <- addLayer(mo_2_2020_02_03, cire_2_2020_03_04 ) 
names(mo_2_2020_02_03) <- paste0(names_S2, "_2mo_2020_02_03")

mo_2_2020_04_05 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2020_04_05.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2020_05_06 <- (mo_2_2020_04_05[[4]]/mo_2_2020_04_05[[8]]) - 1
mo_2_2020_04_05 <- addLayer(mo_2_2020_04_05, cire_2_2020_05_06 ) 
names(mo_2_2020_04_05) <- paste0(names_S2, "_2mo_2020_04_05")

mo_2_2020_06_07 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2020_06_07.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2020_07_08 <- (mo_2_2020_06_07[[4]]/mo_2_2020_06_07[[8]]) - 1
mo_2_2020_06_07 <- addLayer(mo_2_2020_06_07, cire_2_2020_07_08 ) 
names(mo_2_2020_06_07) <- paste0(names_S2, "_2mo_2020_06_07")

mo_2_2020_08_09 <- brick(here("data","DEA_round2", paste0(location,"_2mo_2020_08_09.tif"))) %>% dropLayer(12) #layer 12 is clear pixel count
cire_2_2020_09_10 <- (mo_2_2020_08_09[[4]]/mo_2_2020_08_09[[8]]) - 1
mo_2_2020_08_09 <- addLayer(mo_2_2020_08_09, cire_2_2020_09_10 )
names(mo_2_2020_08_09) <- paste0(names_S2, "_2mo_2020_08_09")


mo_2_2019_10_11_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2019_10_11_SAR.tif")))
names(mo_2_2019_10_11_SAR) <- paste0(names_S1, "_2mo_2019_10_11")

mo_2_2020_12_01_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2020_12_01_SAR.tif")))
names(mo_2_2020_12_01_SAR) <- paste0(names_S1, "_2mo_2020_12_01")

mo_2_2020_02_03_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2020_02_03_SAR.tif")))
names(mo_2_2020_02_03_SAR) <- paste0(names_S1, "_2mo_2020_02_03")

mo_2_2020_04_05_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2020_04_05_SAR.tif")))
names(mo_2_2020_04_05_SAR) <- paste0(names_S1, "_2mo_2020_04_05")

mo_2_2020_06_07_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2020_06_07_SAR.tif")))
names(mo_2_2020_06_07_SAR) <- paste0(names_S1, "_2mo_2020_06_07")

mo_2_2020_08_09_SAR <- brick(here("data","DEA_round2",paste0(location,"_2mo_2020_08_09_SAR.tif")))
names(mo_2_2020_08_09_SAR) <- paste0(names_S1, "_2mo_2020_08_09")

##### combine rasters ####

raster_ready_12m <- addLayer(mo_12, mo_12_SAR)
raster_ready_6m <- addLayer(mo_6_Q41,     mo_6_Q23, 
                            mo_6_Q41_SAR, mo_6_Q23_SAR)
raster_ready_3m <- addLayer(mo_3_Q4,     mo_3_Q1,     mo_3_Q2,     mo_3_Q3, 
                            mo_3_Q4_SAR, mo_3_Q1_SAR, mo_3_Q2_SAR, mo_3_Q3_SAR)
raster_ready_2m <- addLayer(mo_2_2019_10_11,     mo_2_2020_12_01,     mo_2_2020_02_03,     mo_2_2020_04_05,     mo_2_2020_06_07,     mo_2_2020_08_09, 
                            mo_2_2019_10_11_SAR, mo_2_2020_12_01_SAR, mo_2_2020_02_03_SAR, mo_2_2020_04_05_SAR, mo_2_2020_06_07_SAR, mo_2_2020_08_09_SAR)


#### Extract pixel values ####
areas <- read_sf(here("data","FieldData","buffered.geojson"))  %>% st_transform(4326)  %>% 
  dplyr::select(-lon, -lat) 

trainingPoly <- st_read(here("data","FieldData", "newTD_20210727.geojson")) %>%
  dplyr::select(-area) %>% mutate(landcover_level2 = case_when(landcover_level2 == "Cropland irrigated 1" ~ "Cropland irrigated",
                                                               landcover_level2 == "Cropland irrigated 2" ~ "Cropland irrigated",
                                                               landcover_level2 == "Plantation" ~ "Cropland irrigated",
                                                               TRUE ~ landcover_level2),
                                  landcover_code = case_when( landcover_code == 0 ~ "1",
                                                              landcover_code == 7 ~ "2",
                                                              TRUE ~landcover_code)) %>%
  mutate(
    code_level1 = as.character(as.numeric(as.factor(landcover_level1))),
    code_level2 = as.character(as.numeric(as.factor(landcover_level2))),
    PolygonID = seq(1, nrow(.),1),
    DEA_id = 1:nrow(.)) %>% 
  filter(!st_is_empty(.)) %>%
  mutate(point_type = ifelse(is.na(point_type), "Manually drawn", point_type)) #%>% st_transform(crs(four_quarters))

trainingPoly <- st_join(trainingPoly, areas) 

table(trainingPoly$landcover_level2, trainingPoly$area) 
table(trainingPoly$landcover_level2, trainingPoly$landcover_code)
table(trainingPoly$landcover_level2, trainingPoly$point_type )


unique_values <- trainingPoly %>% st_drop_geometry() %>% summarise(landcover_level2  = unique(landcover_level2 ),
                                                                   code = unique(code_level2)) %>% arrange(code)
unique_values # useful for legend in Qgis

#### Function extract and run ####
raster_ready <- raster_ready_6m #change this to the specific composite length

extracted_df <- exact_extract(raster_ready, trainingPoly, force_df = TRUE)

extracted_binded <- bind_rows(extracted_df, .id = "PolygonID") %>% 
  dplyr::select(-coverage_fraction) %>% na.omit() 

TD_df <- merge(extracted_binded, trainingPoly %>% st_drop_geometry(), by = "PolygonID")  %>% select(all_of(names(raster_ready)), "code_level2", "PolygonID")

head(TD_df)
table(TD_df$code_level2)

model_function_ffs <- function(composite, composite_lengths, algorithm){
    set.seed(100)
    polys_split <- initial_split(data = TD_df, prop = .8, strata = code_level2) #prop defines the amount of split #I chose not to split based on collection method, the data was just not good enough to do that
    TD_df_training <- training(polys_split)
    head(TD_df_training)
    
    # model training
    predictors <- names(composite)
    response <- "code_level2"
    
    trainDat <- training(polys_split)

    indices <- CreateSpacetimeFolds(trainDat,
                                    spacevar = "PolygonID",
                                    k=3,
                                    class="code_level2")
    trainDat <- trainDat %>% select(-PolygonID)
    trainDat <- mutate(trainDat, code_level2 = as.factor(code_level2))
    no_cores <- detectCores() - 2
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    
    set.seed = 100
    if(algorithm == "knn"){   # knn does not work in parallel mode
    ctrl <- trainControl(method="cv", 
                         index = indices$index,
                         savePredictions = TRUE,
                         allowParallel= F,
                         number = 5, 
                         verboseIter = TRUE)
    }else{
      ctrl <- trainControl(method="cv", 
                           index = indices$index,
                           savePredictions = TRUE,
                           allowParallel= TRUE,
                           number = 5, 
                           verboseIter = TRUE)}
    
    # num of models = 2x(n-1)^2/2
    n <- nlayers(composite)
    print(2*(n-1)^2/2)
    
  
    if(algorithm == "knn"){
      print("knn model")
      model_ffs <- ffs( trainDat[,predictors],
                        trainDat[,response],
                        method=algorithm, 
                        metric="Accuracy",
                        trControl=ctrl,
                        tuneLength = 5,
                        preProcess = c("center", "scale")
      )
    }else if(algorithm == "nnet"){
      print("nnet model")
      model_ffs <- ffs( trainDat[,predictors],
                        trainDat[,response],
                        method=algorithm, 
                        metric="Accuracy",
                        trControl=ctrl,
                        tuneLength = 10,
                        preProcess = c("center", "scale")
      )
    }else if(algorithm == "svmRadial"){ #svmRadial prediction does not work in ffs mode, so I select the variables using ffs, and do a second training thorugh 'train'
      print("svm model")
      model_ffs_varselect <- ffs( trainDat[,predictors],
                        trainDat[,response],
                        method=algorithm, 
                        metric="Accuracy",
                        trControl=ctrl,
                        importance=TRUE,
                        withinSE = TRUE,
                        tuneLength = 5,
                        na.rm = TRUE ,
                        preProcess = c("center", "scale"))
      
      SVM_radial_vars <- c(model_ffs_varselect$selectedvars ,'code_level2') #here I select only the variables from the ffs output
      trainDat2 <- trainDat %>% select(all_of(SVM_radial_vars)) 
      
      model_ffs <- train(
        code_level2 ~ .,
        data = trainDat2,
        method=algorithm, 
        metric="Accuracy",
        trControl=ctrl,
        importance=TRUE,
        withinSE = TRUE,
        tuneLength = 5,
        na.rm = TRUE ,
        preProcess = c("center", "scale"))
    }    else if(algorithm == "rf"){ #no need to centre or scale RF input data, see discussion: https://stackoverflow.com/questions/8961586/do-i-need-to-normalize-or-scale-data-for-randomforest-r-package
      print( "rf model")
      model_ffs <- ffs( trainDat[,predictors],
                        trainDat[,response],
                        method=algorithm, 
                        metric="Accuracy",
                        trControl=ctrl,
                        importance=TRUE,
                        withinSE = TRUE,
                        tuneLength = 5,
                        na.rm = TRUE 
      )}

    model_ffs <- readRDS(here("output", "Maps", "round_two","Models", paste0(algorithm,"_train_", location, "_", composite_lengths, ".rds")))
    
    saveRDS(model_ffs, here("output", "Maps", "round_two","Models", paste0(algorithm,"_ffs_", location, "_", composite_lengths, ".rds")))
    prediction_ffs_rf <- predict(object = raster_ready, model = model_ffs, progress = "text", cores = no_cores)
    writeRaster(prediction_ffs_rf, here("output", "Maps", "round_two",  "Maps", paste0(algorithm, "_train_", location, "_", composite_lengths, ".tif")),overwrite=TRUE)
}

# composite length models - NOTE: load correct data first
model_function_ffs(raster_ready_12m, "12m", "rf")
model_function_ffs(raster_ready_6m, "6m", "rf")
model_function_ffs(raster_ready_3m, "3m", "rf")
model_function_ffs(raster_ready_2m, "2m", "rf")

# algorithm models
model_function_ffs(raster_ready_6m, "6m", "rf")
model_function_ffs(raster_ready_6m, "6m", "knn")
model_function_ffs(raster_ready_6m, "6m", "svmRadial")
model_function_ffs(raster_ready_6m, "6m", "nnet")

accuracy_function <- function(location, composite_lengths,algo){
  
  prediction_map <- raster(here("output", "Maps", "round_two",  "Maps", paste0(algo,"_ffs", location, "_", composite_lengths, ".tif")))
  
  extracted_df <- exact_extract(prediction_map, trainingPoly, force_df = TRUE)
  
  extracted_binded <- bind_rows(extracted_df, .id = "PolygonID") %>% 
    dplyr::select(-coverage_fraction) %>% na.omit() 
  
  validata_df <- merge(extracted_binded, trainingPoly %>% st_drop_geometry(), by = "PolygonID")  %>% select("value", "code_level2", "PolygonID")
  head(validata_df)

  set.seed(100)
  polys_split <- initial_split(data = validata_df, prop = .8, strata = code_level2)
 
  valiDat <- testing(polys_split) %>% mutate(
    value = as.factor(value),
    code_level2 = as.factor(code_level2)
  )
  
  conmat <- confusionMatrix(data=valiDat$value, reference=valiDat$code_level2)
  
  saveRDS(conmat, here("output", "Maps", "round_two","Accuracy", paste0(algo, "_ffs_", location, "_", composite_lengths, ".rds")))
  
}

# composite length models - NOTE: load correct data first
accuracy_function(location, "12m", "rf")
accuracy_function(location, "6m", "rf")
accuracy_function(location, "3m", "rf")
accuracy_function(location, "2m", "rf")

# algorithm models
accuracy_function(location, "6m", "rf")
accuracy_function(location, "6m", "svmRadial")
accuracy_function(location, "6m", "knn")
accuracy_function(location, "6m", "nnet")

