# Bias Validation

## Goal 

## How-To
Bias validation is ran via `bias_validation.py path/to/toml.input` with an optional `-v` verbose flag. Input files should be stored in `inputs/BiasValidation`, and output will be placed in `outputs`. Below details an exhaustive list of what the toml input can do.
```toml
# First step is to specify the parameters of the full runthrough. This will produce the contours which we wish to validate
# run_name can be any human readable name without whitespace. This will create a pippin input file called BV_FR_run_name.yml in inputs/Pippin. Note that if an input file of that name already exists, the script will ask whether you want to overwrite that script, warning you of the changes. This ensures that we aren't producing the same contours multiple times
[[ full_run.run_name ]]
OMEGA_MATTER = 0.3 # Required. Input cosmology Om 
W0_LAMBDA = -1.0 # Required. Input cosmology w0

# You can specify as many different full runthroughs as you want
[[ full_run.other_run ]]
OMEGA_MATTER = 0.5
W0_LAMBDA = -1.2

# Next we specify the parameters of each validation runthrough. This produces the confidence intervals we care about.  
# As before, run_name_1 can be any human readable name without whitespace. This will create a pippin input file called BV_V_run_name_1.yml in inputs/Pippin. Note that if an input file of that name already exists, the script will ask whether you want to overwrite that script, warning you of the changes. This ensures that we aren't producing the same confidence intervals multiple times
[[ validation.run_name_1 ]]
ompri = 0.3 # Optional, defaults to True. This can either be True or a value between 0 and 1. If true, this will ensure the centre of the gaussian Om prior matches the input Om cosmology. If a value is used instead, then the centre of the Om prior will be set to that
dompri = 0.05 # Optional, defaults to 0.05. This is the width of the gaussian Om prior
OMEGA_MATTER = 0.3 # Required. Can be either a single value or a list of values. This determines the input Om cosmology for the validation runs.
W0_LAMBDA = [-3, -2.5, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0] # Required. Can be either a single value or a list of values. This determines the input w0 cosmology for the validation runs.

# You can specify as many different validation runs as you want
[[ validation.run_name_2 ]]
OMEGA_MATTER = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
W0_LAMBDA = -1.0

# The next step is to actually produce the confidence intervals. This will run the confidence_interval.py script on the output of the validation runthroughs
[ confidence_interval ]
plot = True # Optional, defaults to False. Whether or not to plot the confidence ladder / region
w0 = -1 # Required. Can be either a single value or a list. Specifies for which value of w0 we wish to validate 
Om = 0.3 # Required. Can be either a single value or a list. Specified for which value of Om we wish to validate

# Finally produce the analysis plots. This currently doesn't have any options but is included so that we could add options in the future
[ analyse ]
```

## File Structure

## Inner workings
