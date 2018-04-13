function move_MSC_files_to_surface(subject)
% Script based on Tim's function: post_fc_processing_batch_MSC.m
% CG
% if you don't want to smooth, set smoothnum =0?

% parameters
TR = 2.2;
smoothnum = 2.55;
basedir = '/net/nil-bluearc/GMT/Laumann/MSC/';

% directory information
disp(['Subject #' subject])
bolddir = [basedir '/' subject '/Functionals'];
maskdir = [basedir '/' subject '/subcortical_mask_native_freesurf'];
freesurfdir = ['/net/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/' subject 'edits'];
surfdir = '/net/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR';

%%% FILES TO CHANGE: INPUTS - make these inputs to this file eventually
inputdir = ['/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/MSC05/motor_vol/ttests/']; %input file directory (? I think)
filelist = textread([inputdir 'tstat_files.txt'],'%s');
%%%

% Sample volumes to surface, downsample, and smooth
surffuncdir = [inputdir 'surf_timecourses_native_freesurf/'];
out_filelist = sample_vol_to_surf_func_MSC_cg(filelist,inputdir,smoothnum,surffuncdir,surfdir,subject);

% Smooth data in volume within mask
fcprocess_smooth_volume_wROI_func_cg(filelist,inputdir,maskdir,smoothnum)

% Assumes this is already done:
% Create medial wall mask based on surface projection
% create_medialmask_from_proj_func(tmasklist,surffuncdir,maskdir,smoothnum)

% Create cifti timeseries, normalwall
outputdir = [inputdir '/cifti_timeseries_normalwall_native_freesurf'];
create_cifti_timeseries_manysubs_func_cg(filelist,inputdir,surffuncdir,outputdir,TR,'normalwall',maskdir,out_filelist,smoothnum)

% Create cifti timeseries, smallwall - Only needed for parcellation
%outputdir = [inputdir '/cifti_timeseries_smallwall_native_freesurf'];
%create_cifti_timeseries_manysubs_func(tmasklist,inputdir,surffuncdir,outputdir,TR,'smallwall',maskdir)

end

function out_filelist = sample_vol_to_surf_func_MSC_cg(filelist,funcdir,smoothnum,outputdir,surfdir,subject)
% Function for sampling volumes to the surface following 4 steps:
% 1. Volume-to-surface mapping using ribbon-constrained sampling to native
% surface
% 2. Dilation of data on surface
% 3. Downsample to 32k surface
% 4. Smooth data along surface (if smooth=1)
%
% If a subject is specified as the last variable then tmasklist is treated
% as different sessions of a single subject, rather than different subjects
% TOL, 09/14
%outputdir = [funcdir '/surf_timecourses'];
% CG - 04/16: edited to make more flexible for different input files

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
HEMS = {'L';'R'};
system(['mkdir ' outputdir])

%matlabpool open 8
%parfor s = 1:length(filelist)
for s = 1:length(filelist)
    
    disp(['processing file #' num2str(s) ])
    %subfunc = [funcdir '/' session '/' session '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    subfunc = [funcdir filelist{s}(1:end-9)]; %don't include 4dfp end
    %submask = [bolddir '/' session '/bold' bolduse{1} '/goodvoxels_indiv/' session '_goodvoxels.nii.gz'];
    submask = ['/data/nil-bluearc/GMT/Caterina/goodvoxels_union/' subject '_goodvoxels_union.nii.gz'];
    system(['niftigz_4dfp -n ' subfunc ' ' subfunc]);
    
    for hem = 1:2
        
        midsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        
        surfname = [filelist{s}(1:end-9) '_' HEMS{hem}]; %[filelist{s} '_' HEMS{hem}];
        disp('Map volume to surface')
        system([ workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        disp('Dilate surface timecourse')
        system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
        disp('Deform timecourse to 32k fs_LR')
        system([ workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        disp('Smooth surface timecourse')
        if smoothnum == 0
            system(['cp ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii ' outputdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii']);
        else
            system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii ' num2str(smoothnum) ' ' outputdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'])
        end
        
        system(['rm ' outputdir '/' surfname '.func.gii']);
        system(['rm ' outputdir '/' surfname '_dil10.func.gii']);
        system(['rm ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
        out_filelist{hem,s} =  [surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'];
 
    end
end
%matlabpool close
end

function fcprocess_smooth_volume_wROI_func_cg(filelist,funcdir,maskdir,smoothnum)
% Function volumetrically smooths data within a mask, TOL 09/14
% smooth = 2.55;
% CG - edited it to work for my functions

%[sessions tmasks] = textread(tmasklist,'%s%s');
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

%matlabpool open 6
%parfor s = 1:length(filelist)
for s = 1:length(filelist)

    %sesdir = [funcdir '/' sessions{s}];
    %funcvol = [sesdir '/' sessions{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    funcvol = [funcdir filelist{s}(1:end-9)];
    %system(['niftigz_4dfp -n ' funcvol ' ' funcvol]) - should already be done
    disp(['Volume smooothing, processing file #' num2str(s)])
    system([workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol '_wROI255.nii.gz -roi ' maskdir '/subcortical_mask_LR_333.nii'])

end
%matlabpool close
end

function create_cifti_timeseries_manysubs_func_cg(filelist,funcdir,surffuncdir,outputdir,TR,medialwall,maskdir,out_filelist,smooth_num)
% Create cifti timeseries, TOL 09/14
% Combines masked volume data with surface sampled data to create a cifti
% format file 
% Requires:
% funcvol = nifti format volume of interest
% timename_L = left surface .func.gii file sampled from volume of interest
% timename_L = left surface .func.gii file sampled from volume of interest
% outputdir = output directory
% TR = TR of data for volume timeseries
% medialwall = 'normalwall' or 'smallwall' depending on desired medial wall
% maskdir = directory with medial wall gifti masks

% CG: modified to work for my purposes

%smooth = 2.55;
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
%[subjects tmasks] = textread(tmasklist,'%s%s');
system(['mkdir ' outputdir])

switch medialwall
    case 'smallwall'
        left_mask = [maskdir '/L.atlasroi_noproj.func.gii'];
        right_mask = [maskdir '/R.atlasroi_noproj.func.gii'];
    case 'normalwall'
        left_mask = ['/data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii'];
        right_mask = ['/data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'];
end

matlabpool open 4
parfor s = 1:length(filelist)
%for s = 1:length(filelist) 

    %sub = subjects{s};
    thisFile = filelist{s};
    
    %subdir = [funcdir '/' sub];
    %subfuncvol = [subdir '/' sub '_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255']; 
    subfuncvol = [funcdir '/' filelist{s}(1:end-9) '_wROI255'];
    
    cd(outputdir)
    %Create Cifti timeseries
    disp(['Creating CIFTI timeseries for file number: ' num2str(s)])
    %timename_L = [surffuncdir '/' sub '_L_dil10_32k_fs_LR_smooth' num2str(smooth)];
    %timename_R = [surffuncdir '/' sub '_R_dil10_32k_fs_LR_smooth' num2str(smooth)];
    timename_L = [surffuncdir out_filelist{1,s}(1:end-9)];
    timename_R = [surffuncdir out_filelist{2,s}(1:end-9)];
    
    system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_L '.func.gii'])
    system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_R '.func.gii'])
    %system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' outputdir '/' sub '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smooth) '.dtseries.nii -volume ' subfuncvol '.nii.gz ' maskdir '/subcortical_mask_LR_333.nii -left-metric ' timename_L '.func.gii -roi-left ' left_mask ' -right-metric ' timename_R '.func.gii -roi-right ' right_mask ' -timestep ' num2str(TR) ' -timestart 0'])
    system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' outputdir '/' filelist{s}(1:end-9) '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smooth_num) '.dtseries.nii -volume ' subfuncvol '.nii.gz ' maskdir '/subcortical_mask_LR_333.nii -left-metric ' timename_L '.func.gii -roi-left ' left_mask ' -right-metric ' timename_R '.func.gii -roi-right ' right_mask ' -timestep ' num2str(TR) ' -timestart 0'])
    
end
matlabpool close
end