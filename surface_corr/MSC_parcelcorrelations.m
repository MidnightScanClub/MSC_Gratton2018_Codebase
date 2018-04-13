topdir = '/data/nil-bluearc/GMT/Caterina/';
task = 'mem'
%subjects = 1:10;
subjects = 1:10;

%outputdir_top = [topdir '/TaskFC/FC_Parcels/'];
outputdir_top = [topdir '/TaskFC/FC_Parcels_preGLM/']; %change this and fcdir below as needed

for s = 1:length(subjects)
    
    subject = sprintf('MSC%02d',subjects(s));
    
    tdir_init = [topdir 'TaskFC/FCProc_' subject '_' task '_pass1'];
    %tmasklist = [tdir_init '/COHORTSELECT/NEW_CUT_TMASKLIST.txt'];
    tmasklist = [tdir_init '/COHORTSELECT/NEW_TMASKLIST.txt']; %don't use the "CUT" version; fill empty sessions with nan
    %fcdir = [topdir 'TaskFC/FCProc_' subject '_' task '_pass2'];
    fcdir = [topdir 'TaskFC/FCProc_' subject '_' task '_preGLM_pass2'];

    %Make parcel timecourses
    timecoursedir = [fcdir '/cifti_timeseries_normalwall_native_freesurf'];
    timestem = 'LR_surf_subcort_333_32k_fsLR_smooth2.55';
    watershed_LR = '/data/cn4/laumannt/Parcellation/Parcels_LR.dtseries.nii';
    outputdir = [outputdir_top task]; %[topdir '/TaskFC/FC_Parcels/' task]; 
    [parcel_time, parcel_corrmat parcel_time_concat parcel_corrmat_concat tmask_all] = make_watershed_parcel_timecourse_cifti_func_TASK(tmasklist,timecoursedir,timestem,watershed_LR,outputdir);
    save([outputdir '/' subject '_parcel_timecourse.mat'],'parcel_time','tmask_all')
    save([outputdir '/' subject '_parcel_timecourse_concat_masked.mat'],'parcel_time_concat')
    save([outputdir '/' subject '_parcel_corrmat.mat'],'parcel_corrmat')
    
    %Display correlation matrix
    parcel_correlmat_figmaker_cg(nanmean(parcel_corrmat,3),['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_persess.png']);
    %Display correlation matrix
    parcel_correlmat_figmaker_cg(parcel_corrmat_concat,['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_concat.png']);
    close('all');
end



%%%% And do the same thing for the rest
topdir_rest = '/data/nil-bluearc/GMT/Laumann/MSC/';
task = 'rest'
for s = 1:length(subjects)
    
    subject = sprintf('MSC%02d',subjects(s));
    
    tdir_init = [topdir_rest subject '/Functionals/'];
    tmasklist = [tdir_init 'FCPROCESS_SCRUBBED_UWRPMEAN/corrfile.txt'];
    fcdir = [tdir_init 'FCPROCESS_SCRUBBED_UWRPMEAN'];
    
    %Make parcel timecourses
    timecoursedir = [fcdir '/cifti_timeseries_normalwall_native_freesurf'];
    timestem = 'LR_surf_subcort_333_32k_fsLR_smooth2.55';
    watershed_LR = '/data/cn4/laumannt/Parcellation/Parcels_LR.dtseries.nii';
    outputdir = [outputdir_top task]; %[topdir '/TaskFC/FC_Parcels/' task]; %[dir '/CIMT_MSC02'];
    [parcel_time, parcel_corrmat parcel_time_concat parcel_corrmat_concat tmask_all] = make_watershed_parcel_timecourse_cifti_func_TASK(tmasklist,timecoursedir,timestem,watershed_LR,outputdir,'tmask_from_corrfile');
    save([outputdir '/' subject '_parcel_timecourse.mat'],'parcel_time','tmask_all')
    save([outputdir '/' subject '_parcel_timecourse_concat_masked.mat'],'parcel_time_concat')
    save([outputdir '/' subject '_parcel_corrmat.mat'],'parcel_corrmat')
    
    %Display correlation matrix
    parcel_correlmat_figmaker_cg(mean(parcel_corrmat,3),['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_persess.png']);
    %Display correlation matrix
    parcel_correlmat_figmaker_cg(parcel_corrmat_concat,['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    save_fig(gcf,[outputdir '/' subject '_' task '_parcel_corrmat_concat.png']);
    close('all');
end
