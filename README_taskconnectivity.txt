4.9.2018

Summary of steps for Gratton et al. 2018, Neuron paper

Start with MSC data linked to Gordon et al., 2017 Neuron paper
- Rest data processed through similar steps, as described in that paper
- Below is a description of what was done for task processing, primarily run through MATLAB and FIDL

+++++++++++++++++++++++++
TASK CONNECTIVITY ANALYSES:
+++++++++++++++++++++++++

Steps:

1. Run GLM analysis on data in volume, creating residuals files
- I do this through WUSTL FIDL analysis package
  a. running setupfiles_FIDL_vol.m: this creates residuals script, files needed for FCProcess, links to motion files, etc.
  b. Create GLMs through FIDL: currently this needs to be done by hand - see steps in #2 under TASK ACTIVATION analyses
  c. Create Residual files: run compute_residuals_fidl.csh script
     *** MSC10 mem: needs to be run by hand, since 2nd run is missing

2. Run FC preprocessing, pass 1 [this does first pass of FC process, as described in Power et al., 2014, Neuroimage]

FCPROCESS_MSC_task.m
- EXAMPLE: FCPROCESS_MSC_task('MSC05_motor_DATALIST.txt','/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC05_motor_pass1/','ones',1,'BigBrain264TimOrder_roilist.txt',SUBID,TASK,'Freesurfer/MSC05_FREESURFER.txt','defaults1')

3. Run COHORTSELECT_MSC [this selects good datasets based on FD criteria]
- EXAMPLE: COHORTSELECT_MSC('QC.mat','../MSC05_motor_DATALIST.txt','../MSC05_vcid_edit.txt')
- INPUTS: 0.2 (FD cut-off), 25 (run minimum), 50 (total minimum), 4 (contiguous frame minimum; equates to 5 contiguous frames)
- At this stage, I also update notes on number of frames per subject

4. FCPROCESS_MSC_task.m - PASS 2
- EXAMPLE: FCPROCESS_MSC_task('FCProc_MSC05_motor_pass1/COHORTSELECT/NEW_CUT_DATALIST.txt','/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_MSC05_motor_pass2/','FCProc_MSC05_motor_pass1/COHORTSELECT/NEW_CUT_TMASKLIST.txt',1,'BigBrain264TimOrder_roilist.txt','MSC05','motor','Freesurfer/MSC05_FREESURFER.txt','defaults2')
- NOTE: if running version without task regression, use:
  	FCPROCESS_MSC_task_preGLM2.m
- NOTE2: do to FD high frequency issue (see Fair, in prep), use FDfilt version for:
  	 MSC03, MSC10

5. Run post-FC processing script to put data on the surface, smooth, and calculate correlations
   a. MSC_postFCproc_surface_and_corr.m (use preGLM version for non-residual version)
   b. MSC_parcelcorrelations.m - short version that just remakes correlation matrices (NOTE: this is the one you should run; fills in empty sessions with nans so all are 10 long)

6. Identify frames associated with different task periods:
- in FIDL: export design matrices to a folder called "design_matrices/"
- run: FramesPerCond_MSC(subject,task)
- Update records with frame numbers

7. Run compare_task_rest_mats.m for each task
- Creates session and split-half matrices of each task (and rest each time you do the task for comparison; these are overwritten in successive loops but that should be ok since they should be exactly the same)
- Do this separately for both original (full timeseries) and matched length timeseries
- This also saves out some basic comparisons

8. Run main analyses:
- similarity_figs_SPLITHALF (with both orig and match versions)
- see similar ones for bysess and activation
  for activation analyses, will need to run GLM and create activation maps first (see README_taskactivations.txt), output contrasts of percent signal change, and then extract average activation per parcel using create_task_activation_mats.m
