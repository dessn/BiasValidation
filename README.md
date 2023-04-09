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
[[ reference.run_name ]]
OMEGA_MATTER = 0.3 # Required. Input cosmology Om 
W0_LAMBDA = -1.0 # Required. Input cosmology w0

# Next we specify the parameters of each validation runthrough.
# This produces the confidence intervals we care about.  
# As before, run_name_1 can be any human readable name without whitespace.
# This will create a pippin input file called BV_V_run_name_1.yml in inputs/Pippin.
# Note that if an input file of that name already exists,
# the script will ask whether you want to overwrite that script,
# warning you of the changes.
# This ensures that we aren't producing the same confidence intervals multiple times.
[[ validation.run_name_1 ]]
OMEGA_MATTER = [0.188, 0.38, 0.307, 0.292]
W0_LAMBDA = [-0.783, -1.25, -0.977, -1.02]
# Option to change whether a single biascor is used for every validation
# or whether each validation has its own biascor. Defaults to false
share_biascor = true

# Finally, you can analyse the results
[[ analyse ]]

# Will only grab Pippin output which contains this string.
# Useful if you want to ignore some subset of the validation runs
mask = "SN_OMW_NO"

# Plot the initial contour
# Produces Figure 1
[ analyse.plot_contour ]
# Change the extent of the contour plot, purely aesthetic
extents = [[0.0, 0.6], [-1.6, -0.6]]

# Plot the percentile contour
# This includes the coverage ellipse, and the 150 best-fitting cosmologies
# Produces Figure 2
[ analyse.plot_percentile ]
extents = [[0.0, 0.6], [-1.6, -0.6]]
# Specify every validation run you wish to plot
validation = [0] 

# Plot the process of fitting the coverage ellipse
# This includes the gaussian process fit, and ellipse fit
# Produces Figure 3
[ analyse.plot_GPE ]
validation = [0] 

# Produces Figure 4
# You can't run the same analysis step in the same run through, 
# due to the limitations of TOML
#[ analyse.plot_likelihood ]
#extents = [[0.0, 0.6], [-1.6, -0.6]]
# Plots several validation runs
#validation = [0, 4]
# Optionally plot specific cosmological inputs
#OMEGA_MATTER = 0.138
#W0_LAMBDA = -0.716
# Rename the plot that is produced, defaults to "Percentile.svg"
#name = "Approximate"

# Produces Figure 5
#[ analyse.plot_contour ]
#extents = [[0.0, 0.6], [-1.6, -0.6]]
# Optionally plot specific cosmological inputs
#OMEGA_MATTER = [0.188, 0.38, 0.307, 0.292, 0.145, 0.37, 0.3075, 0.2917]
#W0_LAMBDA = [-0.783, -1.25, -0.977, -1.02, -0.725, -1.204, -0.9765, -1.0205]
#name = "Neyman"

# Plots the coverage ellipse and transformed best fit distributions
# Produces Figure 6
[ analyse.plot_ellipse ]
# Select which validation runs to plot
validation = [0, 1, 2, 3, 4, 5, 6, 7]

# Plots the final cosmologies, which lay on the 68% percentile contour
# Produces Figure 7 and Figure 8
[ analyse.plot_final ]
extents = [[0.0, 0.6], [-1.6, -0.6]]
OMEGA_MATTER = [0.125, 0.37125, 0.308, 0.2918]
W0_LAMBDA = [-0.698, -1.20975, -0.977, -1.02]

```

## Full Runthrough
In order to repeat the analysis in the paper, simply run `bias_validation.py inputs/BiasValidation/lambdacdm.tml`. Warning, if you choose to redo the pippin runthroughs, this will take a significant amount of time (~9 hours on Midway)

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



