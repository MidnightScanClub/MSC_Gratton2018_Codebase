function sample_vol_to_surf_func_MSC_task(tmasklist,funcdir,bolddir,task,smoothnum,outputdir,surfdir,subject)
% function sample_vol_to_surf_func_MSC(tmasklist,funcdir,bolddir,task,smoothnum,outputdir,surfdir,subject)
% Function for sampling volumes to the surface following 4 steps:
% 1. Volume-to-surface mapping using ribbon-constrained sampling to native
% surface
% 2. Dilation of data on surface
% 3. Downsample to 32k surface
% 4. Smooth data along surface
%
% If a subject is specified as the last variable then tmasklist is treated
% as different sessions of a single subject, rather than different subjects
% TOL, 09/14
% CG, edited to work for task data, 01/2017


% some constants
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
HEMS = {'L';'R'};
%smooth = '2.55'; %input now

% 
[sessions tmasks] = textread(tmasklist,'%s%s');
system(['mkdir ' outputdir])

% Account for single subject/many sessions vs. many subjects/single session
if nargin<6
    subjects = sessions;
else
    for s = 1:length(sessions)
        subjects{s} = subject;
    end
end

% Run in parallel - takes a while.
matlabpool open 6
parfor s = 1:length(sessions)
    subject = subjects{s};
    session = sessions{s};
    disp(['processing subject #' num2str(s) ': ' session ])
    
    % functional data for subject
    subfunc = [funcdir '/' session '/' session '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    
    % read out runs to use and mask for good data
    [trash tempruns] = system([ 'awk -F "(" ''{print$2}'' ' funcdir '/' session '/fc.params | awk -F ")" ''{print $1}''' ]);
    bolduse = strread(tempruns,'%s','delimiter',' ');
    submask = [bolddir '/goodvoxels/' session '_b' task bolduse{1} '_goodvoxels.nii.gz']; % CG 
    system(['niftigz_4dfp -n ' subfunc ' ' subfunc]); % convert data to nifti

    % loop through hemispheres
    for hem = 1:2
        
        % set up necessary surfaces
        midsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        
        % run workbench commands for putting volume data on the surface
        surfname = [session '_' HEMS{hem}];
        disp('Map volume to surface')
        system([ workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        disp('Dilate surface timecourse')
        system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
        disp('Deform timecourse to 32k fs_LR')
        system([ workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        disp('Smooth surface timecourse')
        system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii ' num2str(smoothnum) ' ' outputdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'])
         
        % remove files thata are no longer needed
        system(['rm ' outputdir '/' surfname '.func.gii']);
        system(['rm ' outputdir '/' surfname '_dil10.func.gii']);
        system(['rm ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
    end
end
matlabpool close % close pools
