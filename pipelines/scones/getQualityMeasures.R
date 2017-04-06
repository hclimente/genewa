require(tidyverse)
require(magrittr)
require(caret)

getQualityMeasures <- function(pred, truth, modelName){
  c <- confusionMatrix(pred, truth, positive = '1')

  cbind(c$overall %>% as_tibble %>% t,
        c$byClass %>% as_tibble %>% t) %>%
        as_tibble %>%
        set_colnames(. , gsub(" ", "", colnames(.))) %>%
        mutate(model = modelName)

}
