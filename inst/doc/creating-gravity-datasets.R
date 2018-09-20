## ----setup, cache = FALSE, echo = FALSE, message = FALSE, warning = FALSE, tidy = FALSE----
knitr::opts_chunk$set(eval = FALSE)

## ----read----------------------------------------------------------------
#  url <- "http://econ.sciences-po.fr/sites/default/files/file/tmayer/data/col_regfile09.zip"
#  zip <- "col_regfile09.zip"
#  
#  if (!file.exists(zip)) { try(download.file(url, zip)) }
#  try(system("7z e -aos col_regfile09.zip"))
#  
#  library(haven)
#  col_regfile09 <- read_dta("col_regfile09.dta")

## ----isolate-------------------------------------------------------------
#  library(dplyr)
#  data06 <- col_regfile09 %>%
#    filter(year == 2006)

## ----choose--------------------------------------------------------------
#  data06 <- data06 %>%
#    select(iso_o, iso_d, distw, gdp_o, gdp_d, rta, flow, contig, comlang_off, comcur)

## ----complete-cases------------------------------------------------------
#  library(tidyr)
#  gravity_zeros <- data06 %>%
#    drop_na()

## ----scaling-------------------------------------------------------------
#  gravity_zeros <- gravity_zeros %>%
#    mutate(
#      gdp_o = gdp_o / 1000000,
#      gdp_d = gdp_d / 1000000
#    )

## ----no-zeros------------------------------------------------------------
#  gravity_no_zeros <- gravity_zeros %>%
#    filter(flow > 0)

## ----export--------------------------------------------------------------
#  save(gravity_zeros, file = "gravity_zeros.rdata", compress = "xz")
#  save(gravity_no_zeros, file = "gravity_no_zeros.rdata", compress = "xz")

