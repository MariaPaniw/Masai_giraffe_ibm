# Individual-based modeling to assess extinction risks and mitigation for a megaherbivore, the giraffe, in a human-influenced landscape under climate change

This readme.txt file was generated on 2023-06-29 by Maria Paniw


GENERAL INFORMATION

1. Title of Dataset/Repisitory: Data and Analyses from: Extinction risks and mitigation for a megaherbivore, the giraffe, in a human-influenced landscape under climate change

2. Author Information
	A. Principal Investigator Contact Information
		Name: Monica L Bond
		Institution: Doñana Biological Station (EBD-CSIC); Wild Nature Insitute; University of Zurich 
		Address: Seville, 41001 Spain
		Email: monibond@gmail.com

	B. Associate or Co-investigator Contact Information
		Name: Derek E Lee
		Institution:  Wild Nature Insitute; Pennsylvania State University
		Address: Concord, NH, 03301, USA
		Email: derek@wildnatureinstitute.org 

  C. Associate or Co-investigator Contact Information
		Name: Maria Paniw
		Institution: Doñana Biological Station (EBD-CSIC)
		Address: Seville, 41001 Spain
		Email: m.paniw@gmail.com


4. Date of data collection: early 2013 until 2020 for individual giraffe demography, which is the base for the IBM. Rainfall values from 2001 to 2023  

5. Geographic location of data collection: Tarangire Ecosystem in Tanzania (latitude 3.27–4.08°S and longitude 35.73–36.23°E)

6. Information about funding sources that supported the collection of the data: 

Funding for field data collection for the Masai Giraffe Project was provided by Tierpark Berlin, Sacramento Zoo, Tulsa Zoo, Columbus Zoo and Aquarium, The Living Desert Zoo and Gardens, Cincinnati Zoo and Botanical Gardens, Zoo Miami, Toronto Zoo, Como Park Zoo and Conservatory, Roger Williams Park Zoo, Save the Giraffes, and GreaterGood.org.



SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: Please contact authors (A/B in point 2) if you want to use the raw data or code provided here amd to obtain the relevant biological information and prevent misuse of data & code.

2. Links to publications that cite or use the data: NA
   
4. Links to other publicly accessible locations of the data: Paniw, Maria; Bond, Monica (2023). Output of IBM simulations of giraffe populaiton from Tarangire. figshare. Dataset. https://doi.org/10.6084/m9.figshare.23587563.v1

5. Links/relationships to ancillary data sets: Paniw, Maria; Bond, Monica (2023). Output of IBM simulations of giraffe populaiton from Tarangire. figshare. Dataset. https://doi.org/10.6084/m9.figshare.23587563.v1

6. Was data derived from another source? Yes

	Sources (for rainfall data): USGS Earth Resources Observation and Science (EROS) Center (https://earlywarning.usgs.gov/fews) 

7. Recommended citation for this dataset: 
Bond, M. L., Lee, D. E., Paniw, M. (2023). Data from: Extinction risks and mitigation for a megaherbivore, the giraffe, in a human-influenced landscape under climate change. Zenodo. xxxx

DATa & FILE OVERVIEW & WORKFLOW OF R SCRIPTS: 

Scripts are meant to be run in the following order: 

1.	Run core IBM using mean parameters for demographic models: **giraffe_IBM_core_model.R**.

    - a.	*Input*: InitPopGiraffe.csv; rainfall_now.csv
    - b.	*Output*: abund.giraffe.base.csv
 
2.	Plot outputs from core IBM (abundances of different life-cyle stages in 9 communities): **plot_core_abund.R**

    - a.	*Input*: abund.giraffe.base.csv (either generated by user in step 1; or uploaded from https://doi.org/10.6084/m9.figshare.23587563.v1)
    - b.	*Output*: Fig S6 in Supporting Materials
 
3.	Run core IBM including uncertainties in parameters used to construct demographic models: **giraffe_IBM_core_model_PU.R**

    - a.	*Input*: InitPopGiraffe.csv; rainfall_now.csv
    - b.	*Output*: pu.abund.giraffe.base.csv
  
4.	Validate output of the core IBM: **validation.R**

    - a.	*Input*: validation.csv
    - b.	*Output*: Fig. 3 in main text
 
5.	Run sensitivities of changes in abundance to different demographic rates: **giraffe_IBM_core_sensitivity.R**

    - a.	*Input*: InitPopGiraffe.csv; rainfall_now.csv
    - b.	*Output*: sens.giraffe.base.csv

6.	Plot outputs from sensitivity analyses: **plot_sens.R**

    - a.	*Input*: sens.giraffe.base.csv (either generated by user in step 5; or uploaded from https://doi.org/10.6084/m9.figshare.23587563.v1)
    - b.	*Output*: Fig S4 in Supporting Materials

7.	Run IBM perturbing mean parameters for demographic models based on a set of single and combined scenarios of environmental change: **giraffe_IBM_scenarios.R**.

    - a.	*Input*: InitPopGiraffe.csv; rainfall_now.csv; rainfall_10.csv; rainfall_25.csv
    - b.	*Output*: tot.end.scen.csv; ab.time.scen.csv; mean.size.scen.csv; CV.size.scen.csv; ext.df.scen.csv
 
8.	Plot outputs from IBM analysesusing scenarios: **plot_scen.R**

    - a.	*Input*: tot.end.scen.csv; ab.time.scen.csv; mean.size.scen.csv; CV.size.scen.csv; ext.df.scen.csv (either generated by user in step 7; or uploaded from https://doi.org/10.6084/m9.figshare.23587563.v1)
    - b.	*Output*: Figs 4&5 in main text; Figs. S5; S7; S8 in Supporting Materials

9.	Run scenarios IBM including uncertainties in parameters used to construct demographic models (for a subset of main scenarios): **giraffe_IBM_scenarios_PU.R**

    - a.	*Input*: InitPopGiraffe.csv; rainfall_now.csv; rainfall_10.csv; rainfall_25.csv
    - b.	*Output*: tot.end.scen.pu.csv; ab.time.scen.pu.csv; mean.size.scen.pu.csv; CV.size.scen.pu.csv; ext.df.scen.pu.csv

10.	Analyze the contribution of parameter uncertainty to changes in abundances: **giraffe_IBM_scenarios_PU.R**

    - a.	*Input*: tot.end.scen.pu.csv; ab.time.scen.pu.csv; pu.abund.giraffe.base.csv (either generated by user; or uploaded from https://doi.org/10.6084/m9.figshare.23587563.v1)
    - b.	*Output*: Table S4 in Supporting Materials


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 

A detailed description of the methods to develop the IBM is provided in the file (Supporting Materials.pdf). This file will be made availble after the publications of the study. Please contact the authors in the meantime.  


2. Methods for processing the data: 

All data were processed in R

3. Instrument- or software-specific information needed to interpret the data: 

R statistical software, version 4.2.0 (and packages as described in the R scripts)

4. Standards and calibration information, if appropriate: 

5. Environmental/experimental conditions: NA

6. Describe any quality-assurance procedures performed on the data: All R scripts are fully commented and have been checked; all raw data was quality-checked at the study site by data managers.

7. People involved with sample collection, processing, analysis and/or submission: NA 

