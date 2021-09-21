# Scritpt Genlisea lobata #
# 21/09/2021 #

## Carregagar as bibliotecas instaladas ##

library(sp)
library(modleR)
library(raster)

## Importando e lendo a planilha no ambiente R. read.csv é uma função para ler a extensão .csv. ## 
## NO argumento "file" coloque o caminho relativo do arquivo .csv , no arquivo "sep" indique qual ##
## o tipo de separado dos campos (o que separa as colunas). ##

sp_input <- read.csv(file = "./dados/ocorrencias/sp_input_Genlisea_lobata.csv", sep = ",")

## Colocando no formato exigido pelo pacote: species name separated by "_" 
sp_input$species <-
  gsub(x = sp_input$species,
       pattern = " ",
       replacement = "_")

## Carregando as variáveis ambientais

lista_arquivos <- list.files("./dados/raster/Bio_30sec/", full.names = T, pattern = ".tif")

vars_stack <-stack(lista_arquivos)

plot(vars_stack)

## Verificando os pontos nas variáveis ##

par(mfrow = c(1, 1), mar = c(2, 2, 3, 1))
for (i in 1:length(sp_input)) {
  plot(!is.na(vars_stack[[1]]),
       legend = FALSE,
       col = c("white", "#00A08A"))
  points(lat ~ lon, data = sp_input, pch = 19)
}

## modleR função 1 ##

setup_sdmdata_1 <- setup_sdmdata(species_name = unique(sp_input[1]), 
                                 occurrences = sp_input,
                                 lon = "lon",
                                 lat = "lat",
                                 predictors = vars_stack,
                                 models_dir = "./resultados",
                                 partition_type = "crossvalidation",
                                 cv_partitions = 3,
                                 cv_n = 1,
                                 seed = 512,
                                 buffer_type = "mean",
                                 png_sdmdata = TRUE,
                                 n_back = 100,
                                 clean_dupl = TRUE,
                                 clean_uni = TRUE,
                                 clean_nas = TRUE,
                                 #geo_filt = FALSE,
                                 #geo_filt_dist = 0,
                                 select_variables = TRUE,
                                 sample_proportion = 0.5,
                                 cutoff = 0.7)

### Rodando a função do_any para um algoritmo (Maxent) ##

sp_maxent <- do_any(species_name = unique(sp_input[1]),
                    algorithm = "maxnet",
                    #proj_data_folder = "./dados/raster/proj"
                    predictors = vars_stack,
                    models_dir = "./resultados",
                    png_partitions = TRUE,
                    write_bin_cut = FALSE,
                    equalize = TRUE,
                    write_rda = TRUE)

### Rodando a função do_amay para mais de um algoritmo (Maxent) ##

many <- do_many(species_name = unique(sp_input[1]),
                predictors = vars_stack,
                models_dir = "./resultados",
                png_partitions = TRUE,
                write_bin_cut = FALSE,
                write_rda = TRUE,
                bioclim = TRUE,
                domain = FALSE,
                glm = TRUE,
                svmk = FALSE,
                svme = FALSE,
                maxent = FALSE,
                maxnet = FALSE,
                rf = FALSE,
                mahal = FALSE,
                brt = FALSE,
                #proj_data_folder = "./dados/raster/proj"
                equalize = TRUE)

## Usando a função final_model para "unir" as partições geradas por algoritmos em do_any e do_many ##

final_model <- final_model(species_name = unique(sp_input[1]),
                           algorithms = NULL, #if null it will take all the algorithms in disk
                           models_dir = "./resultados",
                           which_models = c("raw_mean",
                                            "bin_mean"),
                           consensus_level = 0.5,
                           uncertainty = FALSE,
                           scale_models = TRUE,
                           overwrite = TRUE)

## Usando a função final_model para "unir" as partições geradas por algoritmos em do_any e do_many ##

ens <- ensemble_model(
  species_name = unique(sp_input[1]),
  occurrences = sp_input,
  which_ensemble = c("average"),
  consensus_level = 0.5,
  which_final = "raw_mean",
  models_dir = "./resultados",
  uncertainty = FALSE,
  overwrite = TRUE
) #argument from writeRaster