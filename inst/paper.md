---
title: 'Gravity: Estimation Methods for Gravity Models in R'
authors:
- affiliation: 1
  name: Anna-Lena Wölwer
- affiliation: 1
  name: Jan Pablo Burgard
- affiliation: 2
  name: Joshua Kunst
- affiliation: 2
  name: Mauricio Vargas
  orcid: 0000-0003-1017-7574
date: "23 August 2018"
output: pdf_document
bibliography: REFERENCES.bib
tags:
- R
- econometrics
- gravity model
- international trade
- migration
affiliations:
- index: 1
  name: Trier University
- index: 2
  name: Pontifical Catholic University of Chile
---

# Summary

Gravity models are used to explain bilateral flows related to the sizes of bilateral partners, 
a measure of distance between them and other influences on interaction costs. 
The underlying idea is rather simple. The greater the masses of two bodies and the smaller the 
distance between them, the stronger their attraction. For a state-of-the-art exposition about 
cross-sectional data see @WoelwerBressleinBurgard2018.

The **gravity** package provides a wrapper of different standard estimation methods that can be quite 
difficult to implement in R [@baser]. By considering the descriptions and code of these 
methods, users get a comprehensive and application-oriented access, can see which methods may be 
suitable for certain research questions or data types, and can extend the code available for their 
specific research projects.

The functions included in this package are designed to be consistent with the Stata code 
used in @Head2014. Beyond offering an **rstats** open alternative to gravity model 
estimation in Stata we also provide cross-system compatibility that eases reproducible research 
for both researchers and students.

The current version of this package relies heavily on the ``rlang`` package [@rlang] which provides 
tools to work with the core language features of base R and the ``tidyverse`` package [@tidyverse]. As a 
result we provide fast model fitting computation while keeping correct handling and consistency 
when facing rank-deficient model matrices that base R handles well.

# References
