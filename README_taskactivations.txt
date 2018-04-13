4.9.2018

+++++++++++++++
TASK ACTIVATION
+++++++++++++++

This is the set of steps used for processing task data through FIDL.
This can be done on either the volume or surface data.
Results are very similar; surface data appears to produce slightly higher peak statistics.

Current steps:
1. Run setupfiles_FIDL.m
   This code sets up file structure and relevant text files

2. Run first level GLM analysis step through FIDL (no easy way to automate this)
    - select GLM -> Design Matrix -> Define single trial design and model
    - Yes, talk to you
    - conc list, all
    - voxel by voxel
    - remove 4
    - all default options
    - accept response parameters
    - response shapes
        - MIXED: sustained have assumed response shape (boxcar); others are modeled with 8 TRs
        - MOTOR: all have assumed response shape (Boynton)
        - MEM: all are modeled with 8 TRs
    - change to 8 time points
    - boxcar for sustained (Boynton for blocked designs)
    - change to mixed_cifti folder for output
    - Exit loop (whenever an option)
    - Boynton contrast
    - No to R2 images
    - Return and run
    - (turn talk off if needed)

3. Check
     design matrix

4. Run t-test/ANOVA scripts that are output for second level analysis
