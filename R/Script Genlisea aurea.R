# Scritpt Genlisea aurea #
# 21/09/2021 #

## Carregagar as bibliotecas instaladas ##

library(sp)
library(modleR)
library(raster)

## Importando e lendo a planilha no ambiente R. read.csv é uma função para ler a extensão .csv. ## 
## NO argumento "file" coloque o caminho relativo do arquivo .csv , no arquivo "sep" indique qual ##
## o tipo de separado dos campos (o que separa as colunas). ##

sp_input <- read.csv(file = "./dados/ocorrencias/sp_input_Genlisea_aurea.csv", sep = ",")

## Colocando no formato exigido pelo pacote: species name separated by "_" 
sp_input$species <-
  gsub(x = sp_input$species,
       pattern = " ",
       replacement = "_")

## Carregando as variáveis ambientais

lista_arquivos <- list.files("./dados/raster/AmbData_Brasil/", full.names = T, pattern = ".asc")

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

#Gerar matriz de correlação 'cor' dos dados acima
cor(tab_sdm)

#Verificação da matriz de correlação entre variáveis
panel.cor <- function(x, y, digits = 2, prefix = "", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.2)
}
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}

pairs(tab_sdm, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor)

#Selecionar as variáveis manualmente 'env.sel' - 'c' combinar as variáveis enumeradas 1,2,3..#
#se na pasta tiver 10 variáveis e apenas quer usar algumas, indicar quais por ordem"
env.sel = env[[c(1,2,3,4)]]
#Selecionar variáveis com a função - na porcentagem adequar à quantidade de pixels dos #
#arquivos de variáveis, se por exemplo tiver 23milhões de pixels, fazer uma regra de 3 para #
#obter uma porcentagem adequada para 5000 pixels que é um valor razoável para um mapa# 
env.sel = select_variables(especies[i], env, percent = 0.0002)
env.sel

#Plotagem dos pontos no mapa raster #
#'legend' para adicionar legendas nos plots #
#'pch' é o formato do ponto no mapa, existe uma sequência no Help que são enumerados 
#de 0 a 25, representando símbolos diversos #
raster::plot(!is.na(env.sel[[1]]), legend = F, ext = ext)
points(sp::SpatialPoints(registros[, c(2, 3)]), bg = (factor(registros$sp)), pch = 19)

#Plotagem da espécie no mapa raster #
raster::plot(!is.na(env.sel[[1]]), legend = F, ext = ext)
points(occs, pch = 19, bg = 1)

## modleR função 1 ##
panel.cor <- function(x, y, digits = 2, prefix = "", ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.2)
}
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
                                 n_back = 50,
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
                           which_models = c("raw_mean_th"),
                           mean_th_par = c("spec_sens"),
                           sensitivity = 0.9,
                           consensus_level = 0.5,
                           uncertainty = FALSE,
                           scale_models = TRUE,
                           png_final = TRUE,
                           overwrite = TRUE)

## Usando a função final_model para "unir" as partições geradas por algoritmos em do_any e do_many ##

ens <- ensemble_model(
  species_name = unique(sp_input[1]),
  occurrences = sp_input,
  which_ensemble = c("average"),
  consensus_level = 0.5,
  which_final = "raw_mean", "bin_consensus",
  png_ensemble = TRUE,
  write_png = TRUE,
  scale_models = TRUE,
  models_dir = "./resultados",
  uncertainty = FALSE,
  overwrite = TRUE
)
