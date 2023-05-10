#!/usr/bin/env Rscript

library(dplyr)
library(data.table)

args <- commandArgs(trailing = TRUE)

ref_file=args[1]
counts_file=args[2]

