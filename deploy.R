#!/usr/bin/env Rscript

#
#  pliu 20211031
#  smcilwain 20221104
#  smcilwain 20251107
#
#  deploy Shiny App to data-viz.it.wisc.edu
#  - if deploy the first time, need to connect via RStudio:
#    - runApp -> publish -> add account -> add server -> connect
#  - please keep rsconnect/ folder so that updated app will be deployed under
#    the same appID and url as well as saving deployment time
#
#  Results:
#  - https://data-viz.it.wisc.edu/content/223/
#


devtools::install_github("cisrlab/PP2A.B56.SLiMs@fimo")

rsconnect::deployApp( appDir  = './',
                      appName = 'itc_predict',
                      account = 'sjmcilwain', ## change to your wisc NetID
                      server  = 'connect.doit.wisc.edu' ,
		      quarto  = FALSE,
		      forceUpdate = TRUE
)
