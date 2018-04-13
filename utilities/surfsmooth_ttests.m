function surfsmooth_ttests(subject,smoothnum,surfaceFile)

%sub_nums = [5]; %[1:10]; %1:10;
%tasks = {'mixed'}; %{'motor','mixed','mem'};


disp('Smooth surface contrasts')

%smoothnum = '2.55';
%subject = 'MSC05';

%workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
workbenchdir = '/data/cn5/caterina/workbench/bin_linux64/';
surfdir = '/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR';
midsurf_LR32k_L = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.L.midthickness.32k_fs_LR.surf.gii'];
midsurf_LR32k_R = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.R.midthickness.32k_fs_LR.surf.gii'];

%dataDir = '/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/MSC05/motor_surf/ttests/';
%surfaceFile = [dataDir '0_Tongue-1_L_Hand-2_R_Hand-3_L_Leg-4_R_Leg_zstat'];

command = [workbenchdir 'wb_command -cifti-smoothing ' surfaceFile '.dscalar.nii ',...
    smoothnum ' ' smoothnum ' COLUMN ' surfaceFile '_sm' smoothnum '.dscalar.nii ',...
    '-left-surface ' midsurf_LR32k_L ' -right-surface ' midsurf_LR32k_R];
system(command);

%%%% NEED TO ADD VOLUME SMOOTHING VERSION too? Tho seems like it would be
%%%% encompassed in the cifti-smoothing command. Check.