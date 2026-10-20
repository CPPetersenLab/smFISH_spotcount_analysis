{\rtf1\ansi\ansicpg1252\cocoartf2870
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 This repository has a workflow in R for analysis of single-molecule FISH data from Guan and Petersen, 2026.\
\
Each folder contains raw data and the analysis code to process this data.  The data is a .csv file where each column is an experimental sample, column names refer to sample IDs, and each row contains a value for transcripts per cell.  The transcripts per cell data was generated using the CPPetersenLab/wholemount_smFISH_quantification pipeline.  \
\
The analysis code described here generates dot plots and ridgeline plots to display each dataset.  It models the data using a hurdle negative binomial model in which non-expressing cells are modeled with a logistic function (hurdle) and expressing cells modeled with negative binomial distribution.  This model was intended to capture the biological variability of smFISH detecting wound-induced genes in whole-mount planarians, in which many cells lack any expression but cells activated for expression display a wide range of target mRNA abundance. The code outputs figures overlaying modeled averages and confidence intervals along with the experimental distributions, and also tables with p-values.\
\
The code is written to read the input file from a working directory that the user must specify when the dataset and analysis is downloaded on their own machine. This input directory also serves as the output directory for plots and tables. }