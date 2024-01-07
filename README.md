# PlasticitySyndromes
Data and code for fitting models used to test for correlated plasticities across multiple behaviors in Johansson et al. 2024. 

**Data Files:**  
mate.switch.data.csv -- mate switch behavioral data for single_behavior_models code  
trip.dist.data.csv -- trip distance behavioral data for single_behavior_models code  
trip.fid.data.csv -- site fidelity behavioral data for single_behavior_models code  
lay.date.data.csv -- laying date behavioral data for single_behavior_models code  
LRS.df -- long term reproductive success data for lay_mate_dhglm code   
annualRS.df -- annual reproductive success data for lay_mate_dhglm code  
dist.fid.dhglm.df -- joint distance/fidelity beahvioral data for dist_fid_dhglm code  
lay.dist.dhglm.df -- joint lay date/distance beahvioral data for lay_dist_dhglm code  
lay.mate.dhglm.df -- joint lay date/mate switch beahvioral data for lay_mate_dhglm code  


**Code Files:**  
single_behavior_models.R -- code for fitting single-behavior models  
dist_fid_dhglm.R -- code for fitting distance-fidelity DHGLM and creating figures   
lay_dist_dhglm.R -- code for fitting lay date-distance DHGLM and creating figures  
lay_mate_dhglm.R -- code for fitting lay date-mate switch DHGLM, reproductive success models, and creating figures  
