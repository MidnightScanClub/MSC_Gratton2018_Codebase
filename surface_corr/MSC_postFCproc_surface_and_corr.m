function MSC_postFCproc_surface_and_corr(subject,task)
% function MSC_postFCproc_surface_and_corr(subject,task)
%%% Script to move task FC proc data to surface for MSC
% CG 4/26/2016
% Based on TOL script, QandD_MSC_processs.m
%
%%% INPUTS
% subject = 'MSC01';
% task = 'motor';
%%%

disp(['*** Processing Subject ' subject ', ' task ' ***']);

% Main directory information
topdir = '/data/nil-bluearc/GMT/Caterina/'; % where task processed/FC processed data lives
timdir = '/data/nil-bluearc/GMT/Laumann/MSC'; % where preprocessed MSC data lives
bolddir = [timdir '/' subject '/Functionals']; 
bolddir_res = [topdir 'SurfTask_analysis/' subject '/' task '_vol/']; % where residuals from task processing live

%%% ACTUAL PROCESSING
% CIFTI information
disp('CIFTI creation')
fcdir = [topdir 'TaskFC/FCProc_' subject '_' task '_pass2']; %where FC processed data lives
maskdir = [timdir '/' subject '/subcortical_mask_native_freesurf']; % subcortical mask for subject
freesurfdir = ['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/' subject 'edits'];  %native Freesurfer data for each subject
surfdir = '/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR'; %group Freesurfer directory
tdir_init = [topdir 'TaskFC/FCProc_' subject '_' task '_pass1']; % FC proc pass one data
tmasklist = [tdir_init '/COHORTSELECT/NEW_CUT_TMASKLIST.txt']; % list of sessions per subject [use all available, not just those that pass cut-off]

%Sample volumes to surface, downsample, and smooth
smoothnum = 2.55; % amount to smooth
TR = 2.2; % data TR
surffuncdir = [fcdir '/surf_timecourses_native_freesurf'];
sample_vol_to_surf_func_MSC_task(tmasklist,fcdir,bolddir_res,task,smoothnum,surffuncdir,surfdir,subject)

%Smooth data in volume within mask
fcprocess_smooth_volume_wROI_func(tmasklist,fcdir,maskdir,smoothnum)

%Create cifti timeseries, normalwall
outputdir = [fcdir '/cifti_timeseries_normalwall_native_freesurf'];
create_cifti_timeseries_manysubs_func(tmasklist,fcdir,surffuncdir,outputdir,TR,'normalwall',maskdir,smoothnum)

%%%%% Parcel timecourses and correlation 
disp('Parcel timecourse and correlation')

%Make parcel timecourses
timecoursedir = [fcdir '/cifti_timeseries_normalwall_native_freesurf'];
timestem = 'LR_surf_subcort_333_32k_fsLR_smooth2.55';
watershed_LR = '/data/cn4/laumannt/Parcellation/Parcels_LR.dtseries.nii';
outputdir = [topdir '/TaskFC/FC_Parcels/' task]; 
[parcel_time, parcel_corrmat parcel_time_concat parcel_corrmat_concat tmask_all] = make_watershed_parcel_timecourse_cifti_func_TASK(tmasklist,timecoursedir,timestem,watershed_LR,outputdir);
save([outputdir '/' subject '_parcel_timecourse.mat'],'parcel_time','tmask_all')
save([outputdir '/' subject '_parcel_timecourse_concat_masked.mat'],'parcel_time_concat')
save([outputdir '/' subject '_parcel_corrmat.mat'],'parcel_corrmat')

%Display correlation matrix
atlas_params = atlas_parameters('Parcels','/data/cn5/caterina/Atlases/Evan_parcellation/');
figure_corrmat_network_generic(mean(parcel_corrmat,3),atlas_params,[-0.4 1]);
save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_persess.png']);
figure_corrmat_network_generic(parcel_corrmat_concat,atlas_params,[-0.4 1]);
save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_concat.png']);