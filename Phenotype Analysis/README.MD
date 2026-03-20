## `1. Fluorcam data` folder

Extracts data from Fluorcam files
   
-Edit the sample list (`Sample list.csv`) to match data
  
## `2. Rapid Light Curves` folder

-Fit the Rapid Light Curve data to `P/I` curves

-Run `Extract_RLC_data.r` which extracts the light curve data from Fluorcam and output it to `outputs/fluorcam_RLC_data.csv`

## 3. Fit parameters to thermal curves

1. Run any `.r` file to fit your parameter(s) of interest to the thermal curve. Depending on the parameter, different models are used.
   - Plots and csv files will be output to the `outputs` folder.
