#install.packages("SPEI")
#Instalar versão SPEI 1.7 ou anterior
#install.packages("readr")
#install.packages("gganimate")
#install.packages("tidyverse")
library(SPEI)
library(tidyverse)
library(gganimate)
library(readr)
library(dplyr)

############### Pontos 0 a 99 ###############	
#Ler o diretorio onde estão todos os arquivos .txt com os pontos 
setwd("G:/Mi unidad/Projeto de Pesquisa Unesp/Índices Hidrometeorológicos/Rstudio/Fenômenos Hidrometeorológicos/Puntos99")
#esta função ajuda ler somente aqueles arquivos com extensão .txt
list.files(pattern = ".txt")
#Importar cada tables com extensão .txt em Rstudio para depois poder visualizar as tabless numa lista 
files <- list.files(recursive = TRUE,pattern = ".txt",full.names = TRUE)
tables <- lapply(files,function(x)read.table(x,header=TRUE,sep="\t"))
#Criar um ciclo for que vai recorrer todas as tabless da lista anterior criada
#em cada tables vai ser calculado variavéis do SPEI, SPI e Balançõ Hídrico
for (i in 1:99) 
{
  #Indicar no script a variavél PET para toda a tables. 
  tables[[i]]$PET <- (tables[[i]]$PET)
  #Crear a variavél de Balanço Hídrico (BAL)
  tables[[i]]$BAL <- tables[[i]]$PRCP - tables[[i]]$PET
  #Criar o Índice Padronizado de Precipitação e Evapotranspiração (SPEI) para escalas mensais, trimestrais e semestrais.
  #SPEI escala Mensal
  spei1 <- spei(tables[[i]][, "BAL"], 1)
  spei1 <- as.numeric(spei1$fitted)
  tables[[i]]$SPEIMen <- spei1
  #SPEI escala Trimestral
  spei3 <- spei(tables[[i]][, "BAL"], 3)
  spei3 <- as.numeric(spei3$fitted)
  tables[[i]]$SPEITri <- spei3
  #SPEI escala Semestral
  spei6 <- spei(tables[[i]][, "BAL"], 6)
  spei6 <- as.numeric(spei6$fitted)
  tables[[i]]$SPEISem <- spei6
  #Criar o Índice Padronizado de Precipitação (SPI) para escalas mensais, trimestrais e semestrais.
  #SPI escala Mensal
  spi1 <- spi(tables[[i]][, "PRCP"], 1)
  spi1 <- as.numeric(spi1$fitted)
  tables[[i]]$SPIMen <- spi1
  #SPI escala Trimestral
  spi3 <- spi(tables[[i]][, "PRCP"], 3)
  spi3 <- as.numeric(spi3$fitted)
  tables[[i]]$SPITri <- spi3
  #SPI escala Semestral
  spi6 <- spi(tables[[i]][, "PRCP"], 6)
  spi6 <- as.numeric(spi6$fitted)
  tables[[i]]$SPISem <- spi6
  #Visualizar os resultados SPEI na tables inicial
  #Convertir a tables em formato de series de tempo (ts) para poder plottear as imagens
  tables[[i]] <- ts(tables[[i]][, -c(1, 2)], end = c(2021, 09), frequency = 12)
  remove(spei1)
  remove(spei3)
  remove(spei6)
  remove(spi1)
  remove(spi3)
  remove(spi6)
}
#Exportar a tables em formato texto (.txt) com os valores obtidos de Balanço Hídrico, SPEIMen, SPEITri e SPEISem, SPIMen, SPITri e SPISem na mesma pasta ou área 	de trabalho
#write.table(tables, "Punto99SPEI.txt")
############### Pontos 150 a 152 ###############
#Ler o diretorio onde estão todos os arquivos .txt com os pontos 
setwd("G:/Mi unidad/Projeto de Pesquisa Unesp/Índices Hidrometeorológicos/Rstudio/Fenômenos Hidrometeorológicos/Puntos150")
#esta função ajuda ler somente aqueles arquivos com extensão .txt
list.files(pattern = ".txt")
#Importar cada tabela com extensão .txt em Rstudio para depois poder visualizar as tabelas numa lista 
arquivo <- list.files(recursive = TRUE,pattern = ".txt",full.names = TRUE)
tabela <- lapply(arquivo,function(x)read.table(x,header=TRUE,sep="\t"))
#Criar um ciclo for que vai recorrer todas as tabelas da lista anterior criada
#em cada tabela vai ser calculado variavéis do SPEI, SPI e Balançõ Hídrico
for (i in 1:53) 
{
  #Indicar no script a variavél PET para toda a tabela. 
  tabela[[i]]$PET <- (tabela[[i]]$PET)
  #Criar a variavél de Balanço Hídrico (BAL)
  tabela[[i]]$BAL <- tabela[[i]]$PRCP - tabela[[i]]$PET
  #Criar o Índice Padronizado de Precipitação e Evapotranspiração (SPEI) para escalas mensais, trimestrais e semestrais.
  #SPEI escala Mensal
  spei1 <- spei(tabela[[i]][, "BAL"], 1)
  spei1 <- as.numeric(spei1$fitted)
  tabela[[i]]$SPEIMen <- spei1
  #SPEI escala Trimestral
  spei3 <- spei(tabela[[i]][, "BAL"], 3)
  spei3 <- as.numeric(spei3$fitted)
  tabela[[i]]$SPEITri <- spei3
  #SPEI escala Semestral
  spei6 <- spei(tabela[[i]][, "BAL"], 6)
  spei6 <- as.numeric(spei6$fitted)
  tabela[[i]]$SPEISem <- spei6
  #Criar o Índice Padronizado de Precipitação (SPI) para escalas mensais, trimestrais e semestrais.
  #SPI escala Mensal
  spi1 <- spi(tabela[[i]][, "PRCP"], 1)
  spi1 <- as.numeric(spi1$fitted)
  tabela[[i]]$SPIMen <- spi1
  #SPI escala Trimestral
  spi3 <- spi(tabela[[i]][, "PRCP"], 3)
  spi3 <- as.numeric(spi3$fitted)
  tabela[[i]]$SPITri <- spi3
  #SPI escala Semestral
  spi6 <- spi(tabela[[i]][, "PRCP"], 6)
  spi6 <- as.numeric(spi6$fitted)
  tabela[[i]]$SPISem <- spi6
  #Visualizar os resultados SPEI na tabela inicial
  #Convertir a tabela em formato de series de tempo (ts) para poder plottear as imagens
  tabela[[i]] <- ts(tabela[[i]][, -c(1, 2)], end = c(2021, 09), frequency = 12)
  remove(spei1)
  remove(spei3)
  remove(spei6)
  remove(spi1)
  remove(spi3)
  remove(spi6)
}
#Exportar a tabela em formato texto (.txt) com os valores obtidos de Balanço Hídrico, SPEIMen, SPEITri e SPEISem, SPIMen, SPITri e SPISem na mesma pasta ou área 	de trabalho
#write.table(tabela, "Punto150SPEI.txt")
