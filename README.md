# MCSED

## Background
MCSED models the optical, near-infrared and infrared spectral energy distribution (SED) of galactic systems.  In light of the fact that there are so many such codes already publicly available, we describe the motivation for MCSED and highlight areas in which this code stands out from the crowd.  First of all, galaxies over cosmic time span a wide range of parameters related to their baryonic content including total stellar mass, gas and stellar metallcity, dust mass and distribution, and star formation history.  This large variation for the totality of all galaxies makes it extremely difficult to develope a general enough SED fitting code to work for all systems.  Instead, most codes are taylored to work best for specific populations.  

MCSED targets galaxies at cosmic noon (z ~ 1-3) that are selected via their emission lines either in the rest-frame optical or ultraviolet wavelengths.  These sources are drawn from the 3DHST survey (http://3dhst.research.yale.edu/Home.html) as well as the HETDEX survey (http://hetdex.org/).  Initial SED fitting efforts revealed that these systems are forming stars at a rate of 1-100 solar masses per year and have total stellar masses of 10^7-10^10 ([Hagen et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...817...79H)).  

We use a monte carlo algorithm, emcee (http://dfm.io/emcee/current/), to then infer the posteriors for fitted parameters such as stellar mass.

## Installation
To acquire and install this code, simply move to a directory where you would like it stored and type:

        git clone https://github.com/grzeimann/MCSED.git

A directory called "MCSED" will be created containing all of the necessary files for the program.  This is a python based code and does require a few standard python based packages.  All of the packages required can be found in the Anaconda distribution environment.  To install Anaconda, see:
https://docs.anaconda.com/anaconda/install/

## How to Run MCSED
The primary script is run_mcsed.py, which can be called from the command line with input arguments.  To view the input arguments, simply type:

        python run_mcsed.py -h

And you will find a help menu like this.
  
        -f, --filename: File to be read for galaxy data
                        
        -s, --ssp: SSP Models, default "fsps"
                        
        -z, --metallicity: Metallicity for SSP models, 0.019 is solar
                        
        -i, --isochrone: Isochrone for SSP model, e.g. "padova"
                        
        -t, --test: Test script with mock data
                        
        -tf, --test_field: Test filters will match the given input field, default "cosmos"
        
All of the available options for MCSED are found in [config.py](https://github.com/grzeimann/MCSED/edit/master/config.py).  Here we break down the most important of those, which include four key variables:

        ssp = 'fsps'  # options include: 'fsps'
        isochrone = 'padova'  # options include: 'padova'
        sfh = 'empirical'  # options include: 'double_powerlaw', 'empirical'
        dust_law = 'noll'  # options include: 'noll'

The single stellar population (SSP) model is left flexible in MCSED as it has a strong weight on the inferred parameters.  We also allow the user to define the isochrone applied within a given SSP model.  To build a composite stellar population the user can set the type of star formation history as well as dust attenuation law to go from instrinsic emission to observed flux.  We are still in the early phases of development and so not many options for either of these four key variables are available.  However, as these are flexible and modular in nature, we expect this to have a great impact on any experimental design.

Another key feature in SED fitting is the metallicity of the SSPs.  We offer two options: fixed metallicity or a relationship between stellar mass and stellar metallicty as given by [Ma et al. (2016)](http://adsabs.harvard.edu/abs/2016MNRAS.456.2140M). (NOTE: only fixed metallicity is currently available.

        # Use metallicity-mass relationship from Ma et al. 2016
        metallicity_mass_relationship = False  # set to True for this option
        # Fixed metallicity of SSP models if fit_metallicity is False
        metallicity = 0.0006  # for fixed metallicity


## Authors

* Greg Zeimann, UT Austin, grzeimann@gmail.com or gregz@astro.as.utexas.edu
