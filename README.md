# MCSED
## Authors

* Greg Zeimann, UT Austin, grzeimann@gmail.com or gregz@astro.as.utexas.edu

## Background
MCSED models the spectral energy distribution (SED) of weakly resolved or unresolved galaxies.  We use parametrized star formation histories (SFHs) and dust attenuation laws in combination with single stellar population models to create a composite stellar population for our observed integrated SEDs to then infer the stellar mass of the system as well as the variables related to SFH and dust.

## Installation
To acquire and install this code, simply move to a directory where you would like it stored and type:

git clone https://github.com/grzeimann/MCSED.git

A directory called "MCSED" will be created containing all of the necessary files for the program.  This is a python based code and does require a few standard python based packages.  All of the packages required can be found in the Anaconda distribution environment.  To install Anaconda, see:
https://docs.anaconda.com/anaconda/install/

## Example uses
The primary script is run_mcsed.py, which can be called from the command line with input arguments.  To view the input arguments, simply type:

python run_mcsed.py -h

And you will find a help menu like this.

usage: run_mcsed_fit.py [-h] [-f FILENAME] [-s SSP] [-z METALLICITY]
                        [-i ISOCHRONE] [-t] [-tf TEST_FIELD]

MCSED

optional arguments:
    -h, --help            show this help message and exit
  
  -f FILENAME, --filename FILENAME
  
                        File to be read for galaxy data
                        
  -s SSP, --ssp SSP     SSP Models, default fsps
  
  -z METALLICITY, --metallicity METALLICITY
  
                        Metallicity for SSP models, 0.02 is solar
                        
  -i ISOCHRONE, --isochrone ISOCHRONE
  
                        Isochrone for SSP model, e.g. padova
                        
  -t, --test            Test script with fake data
  
  -tf TEST_FIELD, --test_field TEST_FIELD
  
                        Test filters will match the given field
