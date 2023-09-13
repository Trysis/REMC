# Replica Exchange Monte Carlo algorithm reimplementation
The goal of this project was to reproduce the REMC algorithm: ![Replica Exchange Monte Carlo algorithm](https://pubmed.ncbi.nlm.nih.gov/17875212/)  

This project was made for the M2 Bioinformatics course in Université Paris Cité.  

## Authors

Roude JEAN MARIE  

## Run

We need to be located in the `src` directory and execute the specified commands :    

usage : ```python REMC [fastafile] [-o output_directory]```  
	```python REMC.py --help  # to get all available options```
	

An animation will be generated after the script is executed, and some output files such as .gif files and energy.png for each replica will be generated in the `./out/` directory if `-o` option value was not specified.
A Pymol file will also be created in the output directory.  


## Installation requierement

Clone the repository :  

```git clone https://github.com/Trysis/REMC.git```  

Move the freshly created repository :  

```cd ./REMC/```  

Import in conda environment :  

`conda env create -f environment.yml`  

Then activate environment :  

`conda activate BI`  

## REMC Algorithm

Given a protein sequence the goal is to apply a set of move to obtain an optimal conformation that will minize the energy of the conformation. A temperature is attributed to a set of replicas and are exchanged given a certain probability during the algorithm.  

Only 2D coordinates were used in this project.  

## Results

VSHD, Pull and Hybrid neighbourhood moves set have been implemented.  
Conformational change animation will be generated after the script execution.  

With the VSHD neighbourhood moves set :  
![front-page](out/run_0/r1_len=16.gif)  

Some energy value during the change in conformation :  
![front-page](out/run_0/energy_r0_len=16.png)

Some 2D to 3D graphical representation with Pymol :  
![front-page](out/run_0/test_pymol.png)

