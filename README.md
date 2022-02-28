# Bias Validation

## Goal 

## How-To
Bias validation is ran via `bias_validation.py path/to/toml.input` with an optional `-v` verbose flag. Input files should be stored in `inputs/BiasValidation`, and output will be placed in `outputs`. Below details an exhaustive list of what the toml input can do.
```toml
# First step is to specify the parameters of the full runthrough.
# This will produce the contours which we wish to validate
# run_name can be any human readable name without whitespace.
# This will create a pippin input file called BV_FR_run_name.yml in inputs/Pippin.
# Note that if an input file of that name already exists,
# the script will ask whether you want to overwrite that script,
# warning you of the changes.
# This ensures that we aren't producing the same contours multiple times
[[ full_run.run_name ]]
OMEGA_MATTER = 0.3 # Required. Input cosmology Om 
W0_LAMBDA = -1.0 # Required. Input cosmology w0

# You can specify as many different full runthroughs as you want
[[ full_run.other_run ]]
OMEGA_MATTER = 0.5
W0_LAMBDA = -1.2

# Next we specify the parameters of each validation runthrough.
# This produces the confidence intervals we care about.  
# As before, run_name_1 can be any human readable name without whitespace.
# This will create a pippin input file called BV_V_run_name_1.yml in inputs/Pippin.
# Note that if an input file of that name already exists,
# the script will ask whether you want to overwrite that script,
# warning you of the changes.
# This ensures that we aren't producing the same confidence intervals multiple times.
[[ validation.run_name_1 ]]
# Optional, defaults to true.
# This can either be true or a value between 0 and 1.
# If true, this will ensure the centre of the gaussian Om prior matches the input Om cosmology.
# If a value is used instead, then the centre of the Om prior will be set to that
ompri = 0.3 
dompri = 0.05 # Optional, defaults to 0.05. This is the width of the gaussian Om prior
# Required. Can be either a single value or a list of values.
# This determines the input Om cosmology for the validation runs.
OMEGA_MATTER = 0.3 
# Required. Can be either a single value or a list of values.
# This determines the input w0 cosmology for the validation runs.
W0_LAMBDA = [-3, -2.5, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0] 

# You can specify as many different validation runs as you want
[[ validation.run_name_2 ]]
OMEGA_MATTER = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
W0_LAMBDA = -1.0

# The next step is to actually produce the confidence intervals.
# This will run the confidence_interval.py script on the output of the validation runthroughs
[ confidence_interval ]
plot = true # Optional, defaults to False. Whether or not to plot the confidence ladder / region
# Required. Can be either a single value or a list.
# Specifies for which value of w0 we wish to validate 
w0 = -1 
# Required. Can be either a single value or a list.
# Specified for which value of Om we wish to validate
Om = 0.3 

# Finally produce the analysis plots.
# This currently doesn't have any options but is included so that we could add options in the future
[ analyse ]
```

## File Structure

On midway, you can find the bias validation code at `/project2/rkessler/SURVEYS/DES/USERS/parmstrong/dev/bias_validation`, hereafter called `$BIAS`.

$BIAS

&emsp;bias_validation.py - The main bias validation file, which can be called alongside an input file to perform bias validations.

&emsp;README.md - The file you are currently reading.

&emsp;files/ - Contains all the base files, input files, and scripts needed to run `bias_validation.py`

&emsp;&emsp;base_files/ - Contains SNANA and BBC base files, used in the pippin scripts run by `bias_validation.py`. This runs wfit with ompri = 0.3, dompri = 0.05, wmin - -2.0 and wmax = 0.0

&emsp;&emsp;&emsp;SALT2mu.input - BBC input file used for a full run with Ompri.

&emsp;&emsp;&emsp;SALT2mu_master.input - BBC input file used for  a validation run. Contains `WFITMUDIF_OPT: $WFIT` which get replaced by `bias_validation.py` based on the simulation inputs for that validation run.

&emsp;&emsp;&emsp;SALT2mu_no_ompriu.input - BBC input file used for a full run without Ompri. This runs wfit with ompri = 0.3, dompri = 5, wmin = -2.0 and wmax = 0.0

&emsp;&emsp;&emsp;des_3yr.yml - LCFIT input file used by all pippin runs

&emsp;&emsp;&emsp;sn_ia_salt2_g10_des3yr - SIM input file used by all pippin runs

&emsp;&emsp;input_files/ - Master pippin input files

&emsp;&emsp;&emsp;full_run.yml - Master pippin input file for a full runthrough. Contains $W0_LAMBDA, and $OMEGA_MATETR which get replaced based on the simulation input.

&emsp;&emsp;&emsp;validation.yml - Master pippin input file for a validation run. Contains $DDATA_DIR, $DESSIM, and $BCOR which get replaced based on the simulation input

&emsp;&emsp;scripts/ - Contains scripts ran by `bias_validation.py`

&emsp;&emsp;&emsp;confidence_intervals.py - Given the output of a validation run, calculate the confidence intervals and plot the confidence ladder

&emsp;&emsp;&emsp;contour.py - Given the output of a full runthrough, and the confidence intervals, plot the contour with the confidence intervals.

&emsp;inputs/ - Contains inputs to both `bias_validation.py` and pippin.

&emsp;&emsp;BiasValidation/ - Contains user created bias validations.

&emsp;&emsp;Pippin/ - Contains the pippin input files generated by `bias_validation.py`. Stored here so the can be easily reused to save time.

&emsp;&emsp;&emsp;SALT2mu/ - Contains the BBC input files generated by `bias_validation.py`. Stored in a central location where pippin can find them.

&emsp;outputs/ - Contains the output of `bias_validation.py`


## Inner Workings
This section will describe in as much detail as is feasable, how the `bias_validation.py` script actually works.

### Full Runthrough
The first step is to prepare and run the full runthroughs. This is done by editing the `files/input_files/full_run.yml` pippin input, based on the options chosen in the toml input file. The only options which changer ar the input cosmology (i.e `W0_LAMBDA` and `OMEGA_MATTER`). This pippin file generates 75 different  
