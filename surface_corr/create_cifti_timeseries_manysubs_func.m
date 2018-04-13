function create_cifti_timeseries_manysubs_func(tmasklist,funcdir,surffuncdir,outputdir,TR,medialwall,maskdir,smooth)
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
% CG - edit to add smooth to input


%smooth = 2.55;
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
[subjects tmasks] = textread(tmasklist,'%s%s');
system(['mkdir ' outputdir])

% which medial wall to use
% use normalwall unless doing parcellations
switch medialwall
    case 'smallwall'
        left_mask = [maskdir '/L.atlasroi_noproj.func.gii'];
        right_mask = [maskdir '/R.atlasroi_noproj.func.gii'];
    case 'normalwall'
        left_mask = ['/data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii'];
        right_mask = ['/data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'];
end

% parallel pool, loop across subjects
matlabpool open 4
parfor s = 1:length(subjects)
    
    % subject information
    sub = subjects{s};
    subdir = [funcdir '/' sub];
    subfuncvol = [subdir '/' sub '_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255'];  % CG: this is for task residual files
    cd(outputdir)
    
    %Create Cifti timeseries
    disp(['Creating CIFTI timeseries for subject number: ' num2str(s)])
    timename_L = [surffuncdir '/' sub '_L_dil10_32k_fs_LR_smooth' num2str(smooth)];
    timename_R = [surffuncdir '/' sub '_R_dil10_32k_fs_LR_smooth' num2str(smooth)];
    
    system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_L '.func.gii'])
    system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_R '.func.gii'])
    system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' outputdir '/' sub '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smooth) '.dtseries.nii -volume ' subfuncvol '.nii.gz ' maskdir '/subcortical_mask_LR_333.nii -left-metric ' timename_L '.func.gii -roi-left ' left_mask ' -right-metric ' timename_R '.func.gii -roi-right ' right_mask ' -timestep ' num2str(TR) ' -timestart 0'])
    
end
matlabpool close
