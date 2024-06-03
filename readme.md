# Introduction

This repository contains the code, data, and supplementary materials for the paper "Exploring the Evolution of the Gene Ontology and its Impact on Enrichment Analysis" presented at the 2024 International Conference on Biomedical Ontologies (ICBO). 

## Table of Contents
- [Introduction](#introduction)
- [pre-requistes](#prerequistes)
- [running basic functions](#running basic functions)



## Introduction 
This project has the following aims:
- Extract Gene Ontology(GO) and Gene Ontology Annotation(GOA) files closest to the select timepoint
- Provide information related to the GO and GOKB, such as: terms & their namespaces, common ancestors of a given pair of terms, IC value of a particular GO term using Resnik method etc.  
- Gene Ontology enrichment analysis 
- analyze the influence of GO evolution on GO enrichment analysis and visualize the differences. 

## prerequisites
This project uses:
- python 3.8.17
- requests==2.29.0
- networkx==3.1
- pandas==2.0.3
- numpy==1.24.4
- obonet==1.0.0
- statsmodels==0.14.1
- sklearn==0.0.post7
- tqdm==4.65.0

To ensure all necessary packages have been installed, run `pip install -r requirements.txt`

## running basic functions