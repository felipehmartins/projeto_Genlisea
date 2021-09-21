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

