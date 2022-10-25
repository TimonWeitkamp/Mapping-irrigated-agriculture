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
library(tmap)
library(ggalluvial)
library(cvms)
library(gt)
library(MLeval)
library(pROC)
library(glue)

image_folder <- ".../img/" 

pixel_count <- function(layer){
  table <- tabularaster::as_tibble(layer) %>%
    count(Prediction = cellvalue) %>% 
    mutate(area = n * 100,
           total_n = sum(n)) %>%  na.omit() %>% #not sure why there is NA
    mutate(total_area = sum(area),
           # class_weight = area/total_area,
           Prediction = as.character(Prediction)) %>% 
    mutate(perc = n/total_n*100)
  return(table)
}

# following Olofsson et al 2014: Good practices for estimating area and assessing accuracy of land change (https://doi.org/10.1016/j.rse.2014.02.015)
unbiased_estimates <- function(composite_length,  area, algo){
  
  prediction <- raster(here("output", "Maps", "round_two",  "Maps", paste0(algo,"_ffs", area, "_", composite_length, ".tif")))
  
  total_area <- tabularaster::as_tibble(prediction) %>%
    count(Prediction = cellvalue) %>% 
    mutate(area = n * 100) %>%  na.omit() %>% #not sure why there is NA
    mutate(total_area = sum(area),
           class_weight = area/total_area,
           Prediction = as.character(Prediction))
  
  result <- readRDS(here("output", "Maps", "round_two","Accuracy", paste0(algo, "_ffs_", area, "_", composite_length, ".rds")))

  conf_mat = result$table
  col_sums <- conf_mat %>% as.data.frame() %>% group_by(Reference) %>% summarise(col_sum= sum(Freq))
  row_sums <- conf_mat %>% as.data.frame() %>% group_by(Prediction ) %>% summarise(row_sum= sum(Freq))
  
  matrix_long <- conf_mat %>% as.data.frame()
  
  conf_mat_weighted <- left_join(total_area, row_sums) %>% 
    left_join(., matrix_long) %>%
    mutate(weighted_value = class_weight*(Freq/row_sum)) %>% 
    group_by(Reference) %>% mutate(col_N_weighted= sum(weighted_value)) %>% ungroup() %>%
    group_by(Prediction) %>% mutate(row_N_weighted= sum(weighted_value)) %>% 
    mutate(class_area_estimate = col_N_weighted*total_area,
           perc_of_total = class_weight  *100) %>% 
    mutate(part1 = ((class_weight*weighted_value - weighted_value^2)/(row_sum-1)) ) %>%
    group_by(Reference) %>%
    mutate(part2 = sum(part1),
           st_error_of_area = sqrt(part2),
           st_error_of_area_m2 = st_error_of_area*total_area,
           CI_95 = st_error_of_area_m2*1.96,
           CI_accuracy = st_error_of_area*1.96) 
  
  area_comparison <- conf_mat_weighted %>% select(Reference, Prediction, area, class_area_estimate) %>% filter(Reference == Prediction) %>% 
    mutate(area = round(area * 0.0001,0),
           class_area_estimate = round(class_area_estimate * 0.0001, 0)) %>%
    mutate(difference = area - class_area_estimate,
           diff_precent = round(difference/class_area_estimate *100,1),
           algo = algo) %>% select(-Prediction, -difference)
  

  class_names_8classes = c("Built-up area", "Irrigated agriculture", "Rainfed agriculture", "Dense vegetation", "Grassland", "Light vegetation",  "Water", "Wetland")
  
  area_comparison$Reference <- plyr::mapvalues(area_comparison$Reference, from = c("1", "2", "3", "4", "5", "6","7" ,"8" ), to = class_names_8classes)

  OA_unbiassed <- conf_mat_weighted %>% select(Prediction, Reference, weighted_value ) %>% filter(Prediction == Reference)%>% ungroup() %>% 
    summarise(OA_unbiassed = round(sum(weighted_value)*100,2),
              algo = algo)
  
  unbiassed_UPA <- conf_mat_weighted %>% select(Prediction, Reference, weighted_value, col_N_weighted, row_N_weighted, CI_accuracy  ) %>% filter(Prediction == Reference) %>% ungroup() %>% 
    mutate(UA_unbiassed = round(weighted_value/row_N_weighted*100,2),
           PA_unbiassed = round(weighted_value/col_N_weighted*100,2)) %>%
    select(Reference, UA_unbiassed, PA_unbiassed,CI_accuracy)
  
  unbiased_areas <- conf_mat_weighted %>% select(Reference, class_area_estimate , CI_95) %>% unique() %>%
    mutate(class_area_estimate_ha = round(class_area_estimate * 0.0001,0),
           CI_95_ha = round(CI_95*0.0001,0)) %>% select(Reference, class_area_estimate_ha , CI_95_ha)
  
  summary_unbiased <- left_join(unbiased_areas, unbiassed_UPA) %>% rename(class = Reference) %>% mutate(algo = algo)
  
  summary_unbiased$class <- plyr::mapvalues(summary_unbiased$class, from = c("1", "2", "3", "4", "5", "6","7" ,"8" ), to = class_names_8classes)

  newlist <- list("Area comparison"= area_comparison, 
                  "Unbiased accuracies"= OA_unbiassed, 
                  "Unbiased summary"= summary_unbiased)
  return(newlist)
}

unbiased_accus <-function(area){
  composite_12m <- unbiased_estimates("12m", area)
  composite_6m <- unbiased_estimates("6m", area )
  composite_3m <- unbiased_estimates("3m",area )
  composite_2m <- unbiased_estimates("2m", area )
  
  a <- composite_12m$`Unbiased summary`  %>% #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = "12m", area = area,  composite_12m$`Unbiased accuracies`) 
  b <- composite_6m$`Unbiased summary`  %>% #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = "6m", area = area,  composite_6m$`Unbiased accuracies`) 
  c <- composite_3m$`Unbiased summary` %>% # filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = "3m", area = area,  composite_3m$`Unbiased accuracies`) 
  d <- composite_2m$`Unbiased summary` %>%  #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = "2m", area = area,  composite_2m$`Unbiased accuracies`) 
  e <- rbind(a,b,c,d)
  return(e)
}

unbiased_accus_algo <-function(area, composite){
  composite_12m <- unbiased_estimates(composite, area, "rf")
  composite_6m <- unbiased_estimates(composite, area, "svmRadial" )
  composite_3m <- unbiased_estimates(composite,area, "knn" )
  composite_2m <- unbiased_estimates(composite, area, "nnet" )
  
  a <- composite_12m$`Unbiased summary`  %>% #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = composite, area = area,  composite_12m$`Unbiased accuracies`) 
  b <- composite_6m$`Unbiased summary`  %>% #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = composite, area = area,  composite_6m$`Unbiased accuracies`) 
  c <- composite_3m$`Unbiased summary` %>% # filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>% 
    mutate( composite = composite, area = area,  composite_3m$`Unbiased accuracies`) 
  d <- composite_2m$`Unbiased summary` %>%  #filter(class %in% c("Irrigated agriculture", "Rainfed agriculture")) %>%
    mutate( composite = composite, area = area,  composite_2m$`Unbiased accuracies`)
  e <- rbind(a,b,c)
  return(e)
}
accu_catandica <- unbiased_accus("Catandica")
accu_xai <- unbiased_accus("Xai-Xai")
# accu_chokwe <- unbiased_accus("Chokwe")
# accu_messica <- unbiased_accus("Manica")

accu_catandica <- unbiased_accus_algo("Catandica", "6m")
accu_xai <- unbiased_accus_algo("Xai-Xai", "6m")


oa_accus_plot <- rbind(accu_chokwe , accu_messica )%>% ungroup() %>% select(area, OA_unbiassed ,composite)  %>% unique() %>%
  select(area, composite, OA_unbiassed) %>% rename(OA = OA_unbiassed,
                                                   Composite = composite)
oa_accus_plot <- accu_catandica %>% select(area, OA_unbiassed ,composite, algo)  %>% unique() %>%
  select(area, composite, OA_unbiassed, algo) %>% rename(OA = OA_unbiassed,
                                                   Composite = composite)

oa_accus_plot$Composite <- factor(oa_accus_plot$Composite, levels = c("12m", "6m", "3m", "2m"))
oa_accus_plot$area <- factor(oa_accus_plot$area, levels = c("Manica", "Chokwe"))

ggplot(oa_accus_plot, aes(OA, area, col = algo)) + geom_point() 

tab_oa <- oa_accus_plot %>% gt( groupname_col = "algo") %>%
  tab_options(row_group.as_column = TRUE) %>% 
  tab_footnote(footnote = "Overall accuracy",
               locations = cells_column_labels(
                 columns = OA
               )) 
# tab_oa %>% gtsave(glue("{image_folder}/table_oa.png"), expand = 10)

accus <- rbind(accu_catandica , accu_xai ) %>% #select(area, class,UA_unbiassed, PA_unbiassed ,composite) %>%
  pivot_longer(names_to = "Accuracy", cols = c("UA_unbiassed" ,"PA_unbiassed"))
accus$composite  <- factor(accus$composite , levels = c("12m", "6m", "3m", "2m"))
accus$area <- factor(accus$area, levels = c("Manica", "Chokwe"))

#### Plots and tables ####

accu_catandica %>% group_by(area, OA_unbiassed) %>% 
  mutate(OA_unbiassed = replace(OA_unbiassed, 2:n(), ""))%>% ungroup() %>%
  select(class, composite, UA_unbiassed, PA_unbiassed, OA_unbiassed, algo) %>% 
  # rename(Class = "class",
  #        `Producer Accuracy (%)` = "PA_unbiassed",
  #       `User Accuracy (%)` = "UA_unbiassed",
  #       `Overall accuracy (%)` = "OA_unbiassed") %>%
  rename(Class = "class",
         Producer  = "PA_unbiassed",
         User = "UA_unbiassed",
         Overall  = "OA_unbiassed") %>%
  gt(groupname_col = "algo")%>%
  tab_spanner("Accuracies (%)",3:5)%>%
  tab_options(row_group.as_column = TRUE) %>%
  opt_row_striping(row_striping = FALSE)%>% 
  tab_style(
    style = list(
      cell_fill(color = "lightblue")
    ),
    locations = cells_body(
      # rows = c(2,9,16,23)) #chokwe
      rows = c(2,8,14,20)) #manica
    # 
  ) # %>% gtsave(glue("{image_folder}/table_accu_manica.png"))

accus %>% filter(class == "Irrigated agriculture") %>% ungroup() %>%
  select(area, algo, UA_unbiassed, PA_unbiassed, OA_unbiassed) %>% 
  rename(Area = "area",
         Algorithm = "algo",
         Producer  = "PA_unbiassed",
         User = "UA_unbiassed",
         Overall  = "OA_unbiassed") %>%
  gt(groupname_col = "Area")%>%
  tab_spanner("Accuracies (%)",3:5)%>%
  tab_options(row_group.as_column = TRUE) %>%
  opt_row_striping(row_striping = FALSE)#  %>% gtsave(glue("{image_folder}/table_irri_accu_cat_xai.png"))
 
accu_names <- c("Producer Accuracy", "User Accuracy")
names(accu_names)<- c("PA_unbiassed", "UA_unbiassed")



manica_plot <- ggplot(data = accus %>% filter(area %in% "Manica"))+
  geom_point(aes(x=composite   , y=value  , color= Accuracy) ) +
  geom_point( aes(x = composite, y = OA_unbiassed ),
              shape = 3) +
  coord_flip() +
  facet_wrap(~class) +
  labs(y = "Accuracy value", x = "Composite",
       title =  "User and producer accuracy per class and composite length for Manica",
       caption = "The crosses (+) represent the overall accuracy")+ 
  theme(legend.position = 'bottom') +
  scale_color_hue(labels = c("Producer accuracy", "User accuracy")) +
  theme(legend.title = element_blank())
# ggsave(plot = manica_plot, glue("{image_folder}/manica_first_models.png"), units = "cm", width = 25, height = 17)

chokwe_plot <- ggplot(data = accus %>% filter(area %in% "Chokwe"))+
  geom_point(aes(x=composite   , y=value  , color= Accuracy) ) +
  geom_point( aes(x = composite, y = OA_unbiassed ),
              shape = 3) +
  coord_flip() +
  facet_wrap(~class) +
  labs(y = "Accuracy value", x = "Composite",
       title =  "User and producer accuracy per class and composite length for Chokwe",
       caption = "The crosses (+) represent the overall accuracy")+ 
  theme(legend.position = 'bottom') +
  scale_color_hue(labels = c("Producer accuracy", "User accuracy")) +
  theme(legend.title = element_blank())
# ggsave(plot = chokwe_plot, glue("{image_folder}/chokwe_first_models.png"), units = "cm", width = 25, height = 17)

PA_plot <- ggplot(data = accus %>% filter(Accuracy %in% "PA_unbiassed"))+
  geom_point(aes(x=composite   , y=value  , color= area) ) +
  geom_point( aes(x = composite, y = OA_unbiassed, color = area ),
             shape = 3) +
  coord_flip() +
  facet_wrap(~class) +
  labs(y = "Producer accuracy",
       title =  "Producer accuracy per class, area, composite length, and algorithm",
       caption = "The crosses (+) represent the overall accuracy")+ 
  theme(legend.position = 'bottom')
# ggsave(plot = PA_plot, glue("{image_folder}/pa_first_models.png"), units = "cm", width = 25, height = 17)


UA_plot <- ggplot(data = accus %>% filter(Accuracy %in% "UA_unbiassed"))+
  geom_point(aes(x=composite   , y=value  , color= area) ) +
  geom_point( aes(x = composite, y = OA_unbiassed, color = area ),
              shape = 3) +
  coord_flip() +
  facet_wrap(~class) +
  labs(y = "User accuracy",
       title =  "User accuracy per class, area, composite length, and algorithm",
       caption = "The crosses represent the overall accuracy")+ 
  theme(legend.position = 'bottom')
# ggsave(plot = UA_plot, "C:/Users/timon/OneDrive - Resilience BV/PhD docs/7. Paper 2/Working docs/img/ua_first_models.png", units = "cm", width = 25, height = 17)




