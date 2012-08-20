NNECalc
=======

## NNE Modeling Calculator

### Description

This project is for developing a standard tool for predicting Chlorophyll concentrations, primarily from nitrogen and phosphorous levels. This includes processing queries from the SWAMP database, performing automated QA/QC, aggregating lab data for each field sample, running these data through several regression models, and outputting user friendly data tables and graphics. New models are also being developed and validated with state data using random forest and neural network methods.

More background information on this project can be found [here](http://www.sccwrp.org/ResearchAreas/Nutrients/NutrientCriteriaSupportStudies/BackgroundNutrientNumericEndpoints.aspx).

### Project Station 

* Finishing validating new models (using both the random forest method, and the traditional regression model)
* Preparing to deploy as R package, which will serve as backend to web app

### File Descriptions


* graph_draft.R
Tools to visualize large batches (i.e., data from many different sites) against a prescirbed threshold. Qual2k model graphing tool is functional. Dodds 97/02 model graphing tool still in early stages of development.

* sketch.R
Recorded steps for compiling together data for validating the new California based NNE predictive model. Data are pulled from the SWAMP, SMC, and Prop 50 databases. 

* NNE model.R
Functional. Sketches of scripts involved in modeling chl-a ~ nutrient based on California data.

* NNE_functions.R
All working functions intended to be used in the final product are stored here. This file is slated to be broken down into more managable files.

* test_scripts.R
Test scripts for running preliminary data through NNE_functions.




