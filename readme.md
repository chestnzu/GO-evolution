# Introduction

This repository contains the code, data, and supplementary materials for the paper "Exploring the Evolution of the Gene Ontology and its Impact on Enrichment Analysis" presented at the 2024 International Conference on Biomedical Ontologies (ICBO). 

## Table of Contents
- [Introduction](#introduction)
- [Running docker](#running docker)
- [running basic functions](#running basic functions)



## Introduction 
This project has the following aims:
- Extract Gene Ontology(GO) and Gene Ontology Annotation(GOA) files closest to the select timepoint
- Provide information related to the GO and GOKB, such as: terms & their namespaces, common ancestors of a given pair of terms, IC value of a particular GO term using Resnik method etc.  
- Gene Ontology enrichment analysis 
- analyze the influence of GO evolution on GO enrichment analysis and visualize the differences. 

## Running docker 
- clone this repo to your local
- `docker run -it --entrypoint /bin/bash GOimages`
- `docker run --rm --name gocomparison -v $(pwd)/output:/app/output GOimages <time1> <time2> /app/example_data/covidGroups/`

To ensure all necessary packages have been installed, run `pip install -r requirements.txt`

## running basic functions

## working on
- mount input data set into docker