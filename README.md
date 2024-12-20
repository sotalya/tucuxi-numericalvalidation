# Numerical validation of Tucuxi



## Description
Tucuxi is a software tool designed for clinical pharmacology, particularly for Therapeutic Drug Monitoring (TDM) and Model-Informed Precision Dosing (MIPD). It enables the prediction of drug concentrations based on patient-specific data, including dosage history, covariates, and blood sample measurements. Tucuxi provides a graphical user interface and is built on top of a robust computing engine. The software can be accessed at http://www.tucuxi.ch.

The aim of this code is to numerically validate Tucuxi's predictions by comparing them with NONMEM, a gold standard in popPK software. The following predicted outcomes are used to compare the two softwares: concentrations at sampling time, individual PK parameters and PK profile indicators: AUC<sub>0-24h</sub>, C<sub>min</sub>, C<sub>max</sub>. This code was developed with Python version 3.12.


## Installation 
> git clone https://github.com/sotalya/tucuxi-numericalvalidation.git 

The necessary packages can be installed using the command:
> pip install -r requirements.txt 

A config.ini file is required to use the code. Here is an example: 


>[PATH]
>
>tucucli = /home/tucuxi/tucuxi/tucuxi-core/build/tucucli/tucucli
>
>drugspath_validation = /home/tucuxi/Repository/tucuxi-numericalvalidation/data/Tucuxi_models/
>
>nonmemmodelspath = /home/tucuxi/Repository/tucuxi-numericalvalidation/data/Nonmem_models/
>
>output = /home/tucuxi/Repository/tucuxi-numericalvalidation/output/
>
>queriespath = /home/tucuxi/Repository/tucuxi-numericalvalidation/output/
>
>queryfile = /home/tucuxi/Repository/tucuxi-python-common/dev/templates/query_template.tqf
>
>listtemplate = /home/tucuxi/Repository/tucuxi-python-common/dev/templates/list_template.xml
>
>requesttemplate = /home/tucuxi/Repository/tucuxi-python-common/dev/templates/request_template.xml
>
>inputfilename = /home/tucuxi/Repository/tucuxi-numericalvalidation/data/population_validation_tucuxi_nonmem.csv
>
>
>[BOOLEAN]
>
>runonserver = False
>
>exportprf = False
>
>
>[URL]
>
>computation = http://193.134.218.125:9090/computation


## Running the code
Once the config.init file has been correctly configured and the packages installed, the code can be run directly via the _tucuxi\_numerical\_validation.py_ script. 

## Usage

### Virtual population
The virtual population can be created using the _Population.virtualdrug\_from\_nb\_patients_ function, by specifying the desired number of patients, the number of concentrations and the study design, a dictionary containing information on the amount and number of doses administered and the sampling times for the concentrations. This population can be saved in Excel file using the _Population.virtualdrug\_save\_population\_in\_file_ function. 

```python
pop = Population.virtualdrug_from_nb_patients(nb_patients=400, nb_samples=4,
                                              study_design={"1": {"amount": 10,
                                                                  "nb_doses": 1,
                                                                  "sample_times": [[0.5,1.5], [3,5], [6,10], [18, 30]]},
                                                            "2": {"amount": 30,
                                                                  "nb_doses": 1,
                                                                  "sample_times": [[10, 30]]}})
pop.virtualdrug_save_population_in_file(file_name="virtualpopulation_400patients.csv", folder_name=foldername_init)                                                                  
```

It is also possible to use an Excel file directly using the _Population.virtualdrug\_from\_csv\_file_ function and specifying the run number. An example is given in the /data folder. 

```python
pop = Population.virtualdrug_from_csv_file(configValues['inputfilename'],
                                           model_id="110000")                                                                                                                                
```

### PopPK models
Three files are required to validate a popPK model:
- A runXX_simulations.mod file used to simulate concentrations for virtual patients. It is based on the FORTRAN format of NONMEM. A limit between 1/5 and 5\* the half-life is necessary to generate clinically correct patients. The residual variability is written with an OMEGA to always have the same variability for the same patient and avoid negative concentrations.
- A runXX.mod file, using FORTRAN format and MAXEVAL =0 to estimate the requested parameters.
- A ch.tucuxi.virtualdrug.mdXX.tdd file, written in Tucuxi format for Tucuxi predictions.

PopPK models already validated are available in /data/Nonmem and /data/Tucuxi folders.

It is also possible to create a new popPK model to be validated, by adding the 3 files described above and adding the name of the model in the dictionary in _tucuxi\_numerical\_validation.py_ file at line 65 to indicate the corresponding administration route.


### Run 
To execute the code, use the _run\_numerical\_validation_ function in _tucuxi\_numerical\_validation.py_. An example is given in main on line 1247 of the _tucuxi\_numerical\_validation.py_ file. 

```python
run_numerical_validation(pop, model_id, foldername_init, args)                                                                                                                             
```


### Results
An excel file is created as output, containing information on the patients created, the a priori outcomes and those predicted by NONMEM and Tucuxi, as well as comparison criteria: relative error between the two software or between one software and the apriori values, MPE, RMSE and bioequivalences for AUC<sub>0-24h</sub>, C<sub>min</sub>, C<sub>max</sub>. 

Graphical comparisons of the predictions of Tucuxi versus NONMEM are also created.



## Support
For now on, feel free to contact yann.thoma@heig-vd.ch for support.

## Authors
Anne Ravix, Annie Cathignol, Thierry Buclin, Chantal Csajka, Monia Guidi, Yann Thoma

