function setupfiles_FIDL_vol()
%function setupfiles_FIDL_vol()
% This function sets up relevant files for FIDL and task FC analyses
% 'vol' because this is for running task data in the volume
%
% C Gratton


% Subjects and tasks to analyze
sub_nums = [1:10]; %1:10;
tasks = {'motor','mixed','mem'};

% Paths
topDir = '/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/';
genDir = '/data/nil-bluearc/GMT/';
FCDir = [genDir 'Caterina/TaskFC/'];
dataDir_top = [genDir 'Laumann/MSC/']; %MSC data; similar to what is available through openfMRI

% Loop through subjects
for s = 1:length(sub_nums)
   
    % set up subject information
    sub_id = sprintf('MSC%02d',sub_nums(s));
    subDir = [topDir sub_id '/'];
    if ~exist(subDir)
        mkdir(subDir);
    end
    dataDir = [dataDir_top sub_id '/'];
    disp(['Sub: ' sub_id]);
    
    % Loop through tasks
    for t = 1:length(tasks)
        
        % get location of data and order of VCIDs - note that this varies
        % slightly by task in a couple of cases
        vclist = get_vclist(sub_id,tasks{t});
        
        % set up task and concs dir
        taskDir = [subDir tasks{t} '_vol/'];
        concsDir = [taskDir 'concs/'];
        ttestDir = [taskDir 'ttests/'];
        resDir = [taskDir 'residuals/'];
        motDir = [taskDir 'motion_params/'];
        atlasDir = [taskDir 'atlas/'];
        origfilesDir = [taskDir 'orig_files/'];
        goodvoxDir = [taskDir 'goodvoxels/'];
        if ~exist(taskDir)
            mkdir(taskDir);
            mkdir(concsDir);
            mkdir(ttestDir);
            mkdir(resDir);
            mkdir(motDir);
            mkdir(atlasDir);
            mkdir(origfilesDir);
            mkdir(goodvoxDir);
        end
        
        % get event file information
        eventDir = [genDir 'MSC_glms/' sub_id '/' tasks{t} '/'];
        eventlist = make_event_list(sub_id,tasks{t},vclist);
        
        % write conc list
        conc_fname = [taskDir 'conc_list.list'];
        write_conc_files(conc_fname,concsDir,sub_id,tasks{t},eventDir,eventlist,dataDir,vclist)
        
        % write a glm list
        glm_fname = [taskDir 'glm_list.list'];
        write_glm_list(glm_fname,taskDir,sub_id,tasks{t},vclist)
        
        % write ttest script
        % NOT currently for most coded - t-test analyses done on surface data
        % switch tasks{t}
        %     case 'motor'
        %         write_ttest_script_motor(glm_fname,taskDir,ttestDir,sub_id,vclist);
        %     case 'mixed'
        %         %error()
        %         %write_ttest_script_mixed(glm_fname,taskDir,ttestDir,sub_id); %%% STILL NEED TO EDIT THIS ONE
        %     case 'mem'
        %         %%% STILL NEED TO WRITE THIS ONE
        % end
        
        % write a residuals computation script
        write_residuals_script(glm_fname,conc_fname,tasks{t},taskDir,sub_id,resDir,vclist);
        
        % write other files needed for FCproc
        write_FCproc_scripts(taskDir,resDir,FCDir,tasks{t},sub_id,vclist,motDir,atlasDir,origfilesDir,goodvoxDir)
        
    end
    
    
end
end

function write_glm_list(fname,taskDir,sub_id,task,vclist)

fid = fopen(fname,'w');
fprintf(fid,'number_of_lists:1\n');
fprintf(fid,'number_of_files:10\n');
for i = 1:10
    fprintf(fid,'%s%s_%s_%s.glm\n',taskDir,sub_id,task,vclist{i});
end
fclose(fid);

end

function write_conc_files(fname,concsDir,sub_id,task,eventDir,eventlist,dataDir,vclist)

%%% Write Conc List
fid = fopen(fname,'w');
fprintf(fid,'number_of_lists:2\n');
fprintf(fid,'number_of_files:10\n');
for i = 1:10
    %fprintf(fid,'%s%s_%s%02d.conc\n',concsDir,sub_id,task,i);
    fprintf(fid,'%s%s_%s_%s.conc\n',concsDir,sub_id,task,vclist{i});
end
fprintf(fid,'number_of_files:10\n');
for i = 1:10
    fprintf(fid,'%s%s\n',eventDir,eventlist{i});
end
fclose(fid);


%%% Write Conc files (typical case)
for i=1:10

    % get #s and vcs for bold runs - sometimes exceptions to easy rule
    [bold_runs vclist_runs] = get_boldruns(task,sub_id,i,vclist{i});

    conc_name = sprintf('%s%s_%s_%s.conc',concsDir,sub_id,task,vclist{i});
    fid = fopen(conc_name,'w');
    %fprintf(fid,'number_of_files:2\n');
    %for b = 1:2
    %    fprintf(fid,'%sFunctionals/%s/bold%s%i/%s_b%s%d_faln_dbnd_xr3d_uwrp_atl.4dfp.img\n',...
    %        dataDir,vclist{i},task,b,vclist{i},task,b);        
    %end
    fprintf(fid,'number_of_files:%d\n',length(bold_runs));
    for b = 1:length(bold_runs)
        fprintf(fid,'%sFunctionals/%s/bold%s%s/%s_b%s%s_faln_dbnd_xr3d_uwrp_atl.4dfp.img\n',...
            dataDir,vclist_runs{b},task,bold_runs{b},vclist_runs{b},task,bold_runs{b});        
    end    
    fclose(fid);
end

end

function eventlist = make_event_list(sub_id,task,vclist)

for i = 1:10
    switch task
        case 'mem'
            %eventlist{i} = [sub_id '_' vclist{i} '_session' num2str(i) '_AlphaMemory_eventfile.fidl'];
            eventlist{i} = ['SeparateStimuli/' sub_id '_' vclist{i} '_ses' num2str(i) '_AlphaMemSeparStim_evfile.fidl'];
        case 'motor'
            eventlist{i} = [sub_id '_' vclist{i} '_session' num2str(i) '_motor_eventfile.fidl'];
        case 'mixed'
            eventlist{i} = [sub_id '_' vclist{i} '_IgErr_mixed8tp_session' num2str(i) '_eventfile.fidl'];
    end
end

end

function write_ttest_script_mixed(glmlist_name,procdir,ttestdir,subid)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% make a scratch dir if it does not exist 
scratchDir1 = [procdir 'SCRATCH_ttest1/'];
scratchDir2 = [procdir 'SCRATCH_ttest2/'];
if ~exist(scratchDir1)
    mkdir(scratchDir1);
    mkdir(scratchDir2);
end

% write datfile
write_datfile_mixed(procdir,scratchDir1,scratchDir2,subid); 

% open up file
fname = sprintf([procdir 'compute_ttest_fidl1.csh'],procdir);
fid = fopen(fname,'w');

% begin script
outstr = ['#!/bin/csh\n'...
    'unlimit\n\n'...
    'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
    'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
    'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
    'if(`uname` == Linux) then\n'...
    '\tset dog = `uname -a`\n'...
    '\tset cat = `expr $#dog - 1`\n'...
    '\tif($dog[$cat] == x86_64) then\n'...
    '\t\tset BIN = $BINLINUX64\n'...
    '\telse\n'...
    '\t\tset BIN = $BINLINUX\n'...
    '\tendif\n'...
    'endif\n\n'...
    'if($#argv != 1) then\n\n'];
fprintf(fid,outstr);

for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -scratchdir SCRATCH_ttest1/ -tc  20+21+28+29 20+21+28+29+53+54+61+62 3+4 3+4+36+37 36+37 42 53+54+61+62 9 9+42\n\n',glmlist{i});
end


outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest1.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_Nimage1.4dfp.img" '...
    '-glm_list_file "' procdir 'fidl_ttest1.list" -scratchdir SCRATCH_ttest1/ -var_thresh "1e-10" '...
    '-prepend ' ttestdir ' -glm ' glmlist{i} ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'else\n\n');
fprintf(fid,'nice +19 $BIN/fidl_ttest -driver "%sfidl_ttest1.dat" -clean_up ONLY\n\n',procdir);
fprintf(fid,'endif\n');
fprintf(fid,'#%sglm_list.list',procdir);

fclose(fid);

% open up file
fname = sprintf([procdir 'compute_ttest_fidl2.csh'],procdir);
fid = fopen(fname,'w');

% begin script
outstr = ['#!/bin/csh\n'...
    'unlimit\n\n'...
    'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
    'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
    'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
    'if(`uname` == Linux) then\n'...
    '\tset dog = `uname -a`\n'...
    '\tset cat = `expr $#dog - 1`\n'...
    '\tif($dog[$cat] == x86_64) then\n'...
    '\t\tset BIN = $BINLINUX64\n'...
    '\telse\n'...
    '\t\tset BIN = $BINLINUX\n'...
    '\tendif\n'...
    'endif\n\n'...
    'if($#argv != 1) then\n\n'];
fprintf(fid,outstr);

for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -scratchdir SCRATCH_ttest2/ -tc  42 9\n\n',glmlist{i});
end


outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest2.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_Nimage2.4dfp.img" '...
    '-glm_list_file "' procdir 'fidl_ttest2.list" -scratchdir SCRATCH_ttest2/ -var_thresh "1e-10" '...
    '-prepend ' ttestdir ' -glm ' glmlist{i} ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'else\n\n');
fprintf(fid,'nice +19 $BIN/fidl_ttest -driver "%sfidl_ttest2.dat" -clean_up ONLY\n\n',procdir);
fprintf(fid,'endif\n');
fprintf(fid,'#%sglm_list.list',procdir);

fclose(fid);

end


function write_ttest_script_motor(glmlist_name,procdir,ttestdir,subid,vclist)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% make a scratch dir if it does not exist
scratchDir_paired = [procdir 'SCRATCH_ttest_paired/'];
scratchDir_single = [procdir 'SCRATCH_ttest_single/'];
if ~exist(scratchDir_paired)
    mkdir(scratchDir_paired);
    mkdir(scratchDir_single);
end

% write datfile
write_datfile_motor(procdir,scratchDir_paired,scratchDir_single,subid,vclist);

% open up file
fname = sprintf([procdir 'compute_ttest_paired_fidl.csh'],procdir);
fid = fopen(fname,'w');

% begin script -- PAIRED
outstr = ['#!/bin/csh\n'...
    'unlimit\n\n'...
    'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
    'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
    'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
    'if(`uname` == Linux) then\n'...
    '\tset dog = `uname -a`\n'...
    '\tset cat = `expr $#dog - 1`\n'...
    '\tif($dog[$cat] == x86_64) then\n'...
    '\t\tset BIN = $BINLINUX64\n'...
    '\telse\n'...
    '\t\tset BIN = $BINLINUX\n'...
    '\tendif\n'...
    'endif\n\n'];
fprintf(fid,outstr);

for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -scratchdir SCRATCH_ttest_paired/ -tc  1 2 2+3 2+3+4+5 2+4 3 3+5 4 4+5 5\n\n',glmlist{i});
end

outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest_paired.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_paired_Nimage.4dfp.img" '...
    '-glm_list_file "' procdir 'fidl_ttest_paired.list" -scratchdir SCRATCH_ttest_paired/ -var_thresh "1e-10" '...
    '-prepend ' ttestdir ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sglm_list.list',procdir);

fclose(fid);

%%%% begin script -- SINGLE
fname = sprintf([procdir 'compute_ttest_single_fidl.csh'],procdir);
fid = fopen(fname,'w');

outstr = ['#!/bin/csh\n'...
    'unlimit\n\n'...
    'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
    'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
    'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
    'if(`uname` == Linux) then\n'...
    '\tset dog = `uname -a`\n'...
    '\tset cat = `expr $#dog - 1`\n'...
    '\tif($dog[$cat] == x86_64) then\n'...
    '\t\tset BIN = $BINLINUX64\n'...
    '\telse\n'...
    '\t\tset BIN = $BINLINUX\n'...
    '\tendif\n'...
    'endif\n\n'];
fprintf(fid,outstr);

for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -gauss_smoth 2 -scratchdir SCRATCH_ttest_single/ -tc  1 2 2+3 2+3+4+5 2+4 3 3+5 4 4+5 5\n\n',glmlist{i});
end

outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest_single.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_single_Nimage.4dfp.img" '...
    '-glm_list_file "' procdir 'fidl_ttest_single.list" -scratchdir SCRATCH_ttest_single/ -var_thresh "1e-10" '...
    '-prepend ' ttestdir ' -clean_up\n\n'];
fprintf(fid,outstr);
fclose(fid);



end

function write_FCproc_scripts(procDir,resDir,FCDir,task,subid,vclist,motDir,atlasDir,origfilesDir,goodvoxDir)
%%% some complications with links, etc., to deal with exception cases...


%%% FCParams file
for i = 1:length(vclist)
    [bold_runs,vclist_runs] = get_boldruns(task,subid,i,vclist{i});
    fname = sprintf('%s%s_%s.fcparams',resDir,task,vclist{i});
    fid = fopen(fname,'w');
    str = 'set boldruns = ('; %make flexible for different # of runs
    for b = 1:length(bold_runs)
        str = [str num2str(b) ' ']; % renumber to 1 and 2, since this is how FIDL outputs it
    end
    str = [str ')\n'];
    %fprintf(fid,'set boldruns = (1 2)\n'); % all tasks have 2 runs (?)
    fprintf(fid,str); % all tasks have 2 runs (?)
    fclose(fid);
end

%%% Datalist file
fname = sprintf('%s%s_%s_DATALIST.txt',FCDir,subid,task);
fid = fopen(fname,'w');
for i = 1:10
    fprintf(fid,'%s %s %s%s_%s.fcparams 2.2 0\n',procDir,vclist{i},resDir,task,vclist{i});
end
fclose(fid);

%%% Link in motion parameter files needed for FCproc 
for i = 1:10
    [bold_runs,vclist_runs] = get_boldruns(task,subid,i,vclist{i}); %to make it flexible for exceptions
    for b = 1:length(bold_runs)
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s%s/*.mat %s/%s_b%s%d_faln_dbnd_xr3d.mat',...
            subid,vclist_runs{b},task,bold_runs{b},motDir,vclist{i},task,b);
        system(command);
    end   
    %command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s1/*.mat %s.',subid,vclist{i},task,motDir);
    %system(command);
    %command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s2/*.mat %s.',subid,vclist{i},task,motDir);
    %system(command);
end


%%% Link in anat ave and mpr files for FCproc too
for i = 1:10
    [bold_runs,vclist_runs] = get_boldruns(task,subid,i,vclist{i}); %to make it flexible for exceptions
    for b = 1:length(bold_runs)
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/unwarp_mean/%s_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.img %s/%s_b%d_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.img',...
            subid,vclist_runs{b},vclist_runs{b},atlasDir,vclist{i},b);
        system(command);
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/unwarp_mean/%s_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.ifh %s/%s_b%d_faln_dbnd_xr3d_uwrp_atl_ave.4dfp.ifh',...
            subid,vclist_runs{b},vclist_runs{b},atlasDir,vclist{i},b);
        system(command);
        %command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/atlas/%s_mpr_n1_333_t88.4dfp.img %s/%s_b%s_mpr_n1_333_t88.4dfp.img',...
        %    subid,vclist{i},vclist{i},atlasDir,vclist{i},bold_runs{b});
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/T1/%s_mpr_debias_avgT_333_t88.4dfp.img %s/%s_b%d_mpr_n1_333_t88.4dfp.img',...
            subid,subid,atlasDir,vclist{i},b);
        system(command);
        %command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/atlas/%s_mpr_n1_333_t88.4dfp.ifh %s/%s_b%s_mpr_n1_333_t88.4dfp.ifh',...
        %    subid,vclist{i},vclist{i},atlasDir,vclist{i},bold_runs{b});
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/T1/%s_mpr_debias_avgT_333_t88.4dfp.ifh %s/%s_b%d_mpr_n1_333_t88.4dfp.ifh',...
            subid,subid,atlasDir,vclist{i},b);
        system(command);
    end
    
end

%%% Link in original files for compute defined computations
for i = 1:10
    [bold_runs,vclist_runs] = get_boldruns(task,subid,i,vclist{i}); %to make it flexible for exceptions
    for b = 1:length(bold_runs)
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s%s/%s_b%s%s_faln_dbnd_xr3d_uwrp_atl.4dfp.img %s/%s_b%s%d_faln_dbnd_xr3d_uwrp_atl.4dfp.img',...
            subid,vclist_runs{b},task,bold_runs{b},vclist_runs{b},task,bold_runs{b},origfilesDir,vclist{i},task,b);
        system(command);
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s%s/%s_b%s%s_faln_dbnd_xr3d_uwrp_atl.4dfp.ifh %s/%s_b%s%d_faln_dbnd_xr3d_uwrp_atl.4dfp.ifh',...
            subid,vclist_runs{b},task,bold_runs{b},vclist_runs{b},task,bold_runs{b},origfilesDir,vclist{i},task,b);
        system(command);
    end   
end

%%% Link in goodvoxels files to use
for i = 1:10
    [bold_runs,vclist_runs] = get_boldruns(task,subid,i,vclist{i}); %to make it flexible for exceptions
    for b = 1:length(bold_runs)
        command = sprintf('ln -s /data/nil-bluearc/GMT/Laumann/MSC/%s/Functionals/%s/bold%s%s/goodvoxels_indiv/%s_goodvoxels.nii.gz %s/%s_b%s%d_goodvoxels.nii.gz',...
            subid,vclist_runs{b},task,bold_runs{b},vclist_runs{b},goodvoxDir,vclist{i},task,b);
        system(command);
    end   
end

end

function write_residuals_script(glmlist_name,conclist_name,task,procdir,subid,resDir,vclist)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% also import data from conc list
conc_all = importdata(conclist_name);
conclist = conc_all.textdata(3:12);%%%

% estimates to retain (i.e., the #s for the mean effects)
switch task
    case 'motor'
        retain_estimates = '8 9';
    case 'mixed'
        retain_estimates = '69 70';
    case 'mem'
        retain_estimates = '148 149 150';
end

% open up file
fname = [procdir 'compute_residuals_fidl.csh'];
fid = fopen(fname,'w');

% begin script
outstr = ['#!/bin/csh\n'...
    'unlimit\n\n'...
    'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
    'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
    'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
    'if(`uname` == Linux) then\n'...
    '\tset dog = `uname -a`\n'...
    '\tset cat = `expr $#dog - 1`\n'...
    '\tif($dog[$cat] == x86_64) then\n'...
    '\t\tset BIN = $BINLINUX64\n'...
    '\telse\n'...
    '\t\tset BIN = $BINLINUX\n'...
    '\tendif\n'...
    'endif\n\n'];
fprintf(fid,outstr);

for i = 1:10
    fprintf(fid,'set GLM_FILE = (-glm_file %s)\n',glmlist{i});
    fprintf(fid,'nice +19 $BIN/compute_residuals $GLM_FILE -bold_files %s -retain_estimates %s -root "%s%s_%s_%s_res"\n\n',...
        conclist{i},retain_estimates,resDir,subid,task,vclist{i});
end

fprintf(fid,'# %s',glmlist_name);
fclose(fid);


end



function write_datfile_mixed(procdir,scratchDir1,scratchDir2,subid)


%%%%
% open up file1
fname = sprintf([procdir 'fidl_ttest1.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_Glass_sustained+7_NV_sustained_1\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_1_Glass_sustained+7_NV_sustained_1.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_Glass_sustained_1\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_1_Glass_sustained_1.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 7_NV_sustained_1\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_7_NV_sustained_1.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 0_Glass_startcue+6_NV_startcue_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_0_Glass_startcue+6_NV_startcue_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 0_Glass_startcue_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_0_Glass_startcue_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 6_NV_startcue_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_6_NV_startcue_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 3_Glass_coherent+4_Glass_random+9_Noun+10_Verb_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_3_Glass_coherent+4_Glass_random+9_Noun+10_Verb_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 3_Glass_coherent+4_Glass_random_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_3_Glass_coherent+4_Glass_random_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 9_Noun+10_Verb_3+4\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_9_Noun+10_Verb_3+4.4dfp.img\n',scratchDir1,subid,i);
end

fclose(fid);

%%%%
% open up file
fname = sprintf([procdir 'fidl_ttest2.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'TTEST  := PAIRED_COMPARISON 7_NV_sustained-1_Glass_sustained\n');
for i=1:10
    fprintf(fid,'FIRST  := %s%s_mixed%02d_7_NV_sustained_1.4dfp.img\n',scratchDir2,subid,i);
end
for i=1:10
    fprintf(fid,'SECOND := %s%s_mixed%02d_1_Glass_sustained_1.4dfp.img\n',scratchDir2,subid,i);
end

fclose(fid);

end


function write_datfile_motor(procdir,scratchDir_paired,scratchDir_single,subid,vclist)

%%%% PAIRED VERSION
% open up file
fname = sprintf([procdir 'fidl_ttest_paired.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'TTEST  := PAIRED_COMPARISON 0_Tongue-1_L_Hand-2_R_Hand-3_L_Leg-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_0_Tongue_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 0_Tongue-1_L_Hand-2_R_Hand\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_0_Tongue_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_1_L_Hand+2_R_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 0_Tongue-3_L_Leg-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_0_Tongue_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 1_L_Hand+3_L_Leg-2_R_Hand-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand+3_L_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_2_R_Hand+4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 1_L_Hand-3_L_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_3_L_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 2_R_Hand-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_2_R_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end


fprintf(fid,'TTEST  := PAIRED_COMPARISON 1_L_Hand-2_R_Hand\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_2_R_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fprintf(fid,'TTEST  := PAIRED_COMPARISON 3_L_Leg-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_3_L_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end


fprintf(fid,'TTEST  := PAIRED_COMPARISON 1_L_Hand+2_R_Hand-3_L_Leg-4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand+2_R_Hand_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_motor_%s_3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_paired,subid,vclist{i});
end

fclose(fid);


%%%% SINGLE VERSION
% open up file
fname = sprintf([procdir 'fidl_ttest_single.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 0_Tongue\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_0_Tongue_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_L_Hand\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 2_R_Hand\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_2_R_Hand_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 3_L_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_3_L_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_4_R_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_L_Hand+2_R_Hand\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand+2_R_Hand_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 3_L_Leg+4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_L_Hand+3_L_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand+3_L_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 2_R_Hand+4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_2_R_Hand+4_R_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
end


fclose(fid);


end

function vclist = get_vclist(sub_id,task)
if strcmp(sub_id,'MSC02')
    if strcmp(task,'motor')
        func_file = ['/data/nil-bluearc/GMT/MSC_glms/' sub_id '/' sub_id '_motor_functional_vclist.txt'];
    else
        func_file = ['/data/nil-bluearc/GMT/MSC_glms/' sub_id '/' sub_id '_mixed_memory_functional_vclist.txt'];
    end
else
    func_file = ['/data/nil-bluearc/GMT/MSC_glms/' sub_id '/' sub_id '_functional_vclist.txt'];
end

[vclist_init] = textread(func_file,'%s');

for v = 1:length(vclist_init)
    vclist{v} = ['vc' vclist_init{v}];
end

end

