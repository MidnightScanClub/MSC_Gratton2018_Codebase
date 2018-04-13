function setupfiles_FIDL()
%function setupfiles_FIDL()
%
% Function which sets up files MSC files for FIDL analyses
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
        taskDir = [subDir tasks{t} '_surf/'];
        concsDir = [taskDir 'concs/'];
        ttestDir = [taskDir 'contrasts/'];
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
        switch tasks{t}
            case 'motor'
                write_ttest_script_motor(glm_fname,taskDir,ttestDir,sub_id,vclist);
            case 'mixed'
                write_contrast_script_mixed(glm_fname,taskDir,ttestDir,sub_id,vclist);
            case 'mem'
                if strcmp(sub_id,'MSC10')
                    disp('RUN BY HAND for MSC10, due to missing scene run in session 6');
                else
                    write_contrast_script_mem(glm_fname,taskDir,ttestDir,sub_id,vclist);
                end
        end
        
        % write a residuals computation script
        %write_residuals_script(glm_fname,conc_fname,tasks{t},taskDir,sub_id,resDir,vclist);
        %if strcmp(sub_id,'MSC10')
        %    disp('Remember to edit residuals script by hand to account for missing scene run in session 6');
        %end
        
        % write other files needed for FCproc
        %write_FCproc_scripts(taskDir,resDir,FCDir,tasks{t},sub_id,vclist,motDir,atlasDir,origfilesDir,goodvoxDir)
        
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
        fprintf(fid,'%sFunctionals/%s/bold%s%s_cifti/%s_b%s%s_faln_dbnd_xr3d_uwrp_atl_LR_surf_subcort.dtseries.nii\n',...
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

function write_contrast_script_mixed(glmlist_name,procdir,ttestdir,subid,vclist)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% make scratch directories if they do not exist
scratchDir.paired = [ttestdir 'ttests/SCRATCH_ttest_paired/'];
scratchDir.single = [ttestdir 'ttests/SCRATCH_ttest_single/'];
scratchDir.Fourstim = [ttestdir '4stimXtime/SCRATCH_4stim/'];
scratchDir.endcue = [ttestdir 'endcue/SCRATCH_endcue/'];
scratchDir.startcue = [ttestdir 'startcue/SCRATCH_startcue/'];
scratchDir.Glass = [ttestdir 'Glass/SCRATCH_Glass/'];
scratchDir.NV = [ttestdir 'NV/SCRATCH_NV/'];
scratchDir.tasktime = [ttestdir 'taskXtime/SCRATCH_tasktime/'];
scratchDir.AllCond = [ttestdir 'AllCond/SCRATCH_AllCond/'];
if ~exist(scratchDir.paired)
    mkdir(scratchDir.paired);
    mkdir(scratchDir.single);
    mkdir(scratchDir.Fourstim);
    mkdir(scratchDir.endcue);
    mkdir(scratchDir.startcue);
    mkdir(scratchDir.Glass);
    mkdir(scratchDir.NV);
    mkdir(scratchDir.tasktime);
    mkdir(scratchDir.AllCond);
end

% write datfiles
write_datfile_mixed(ttestdir,scratchDir,subid,vclist); %need vclist?

%%% T-TESTS
% begin script -- PAIRED 1
fname = sprintf('%sttests/compute_ttest_paired_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir %s -tc  42 9\n\n',glmlist{i},scratchDir.paired);
end
outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'ttests/fidl_ttest_paired.dat" '...
    '-output Z_uncorrected  -Nimage_name "fidl_ttest_paired_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir ' scratchDir.paired ' -var_thresh "1e-10" '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sttests/%s',ttestdir,glmlist_name);
fclose(fid);

% begin script -- SINGLE ttest 2
fname = sprintf('%sttests/compute_ttest_single_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir %s -tc  42 9 9+42\n\n',glmlist{i},scratchDir.single);
end
outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'ttests/fidl_ttest_single.dat" '...
    '-output Z_uncorrected  '...
    '-Nimage_name "' ttestdir 'ttests/fidl_ttest_single_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir ' scratchDir.single ' -var_thresh "1e-10" '...
    ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sttests/%s',ttestdir,glmlist_name);
fclose(fid);

% AllCond (vs. Baseline) computation
fname = sprintf('%sAllCond/compute_ttest_single_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir %s -tc  1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29+30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60+61+62+63+64+65+66\n\n',glmlist{i},scratchDir.AllCond);
end
outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'AllCond/fidl_ttest_single.dat" '...
    '-output Z_uncorrected  '...
    '-Nimage_name "' ttestdir 'AllCond/fidl_ttest_single_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir ' scratchDir.AllCond ' -var_thresh "1e-10" '...
    ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sAllCond/%s',ttestdir,glmlist_name);
fclose(fid);


%%%% open up startcue ANOVA file 3
% begin script -- STARTCUE
fname = sprintf('%sstartcue/start_fidl_anova.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 1 2 3 4 5 6 7 8  34 35 36 37 38 39 40 41 -scratchdir %s \n\n',glmlist{i},scratchDir.startcue);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t1,1,1,1,1,1,1,1,1,1 \\\n'...
        '\t2,2,2,2,2,2,2,2,2,2 \\\n'...
        '\t3,3,3,3,3,3,3,3,3,3 \\\n'...
        '\t4,4,4,4,4,4,4,4,4,4 \\\n'...
        '\t5,5,5,5,5,5,5,5,5,5 \\\n'...
        '\t6,6,6,6,6,6,6,6,6,6 \\\n'...
        '\t7,7,7,7,7,7,7,7,7,7 \\\n'...
        '\t8,8,8,8,8,8,8,8,8,8 \\\n'...
        '\t34,34,34,34,34,34,34,34,34,34 \\\n'...
        '\t35,35,35,35,35,35,35,35,35,35 \\\n'...
        '\t36,36,36,36,36,36,36,36,36,36 \\\n'...
        '\t37,37,37,37,37,37,37,37,37,37 \\\n'...
        '\t38,38,38,38,38,38,38,38,38,38 \\\n'...
        '\t39,39,39,39,39,39,39,39,39,39 \\\n'...
        '\t40,40,40,40,40,40,40,40,40,40 \\\n'...
        '\t41,41,41,41,41,41,41,41,41,41 \\\n'...
        '\t)\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'startcue/start_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "start_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.startcue '  -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sstartcue  %s',ttestdir,glmlist_name);
fclose(fid);

%%%% open up endcue ANOVA file 4
% begin script -- ENDCUE
fname = sprintf('%sendcue/end_fidl_anova.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 10 11 12 13 14 15 16 17  43 44 45 46 47 48 49 50 -scratchdir %s \n\n',glmlist{i},scratchDir.endcue);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t10,10,10,10,10,10,10,10,10,10 \\\n'...
        '\t11,11,11,11,11,11,11,11,11,11 \\\n'...
        '\t12,12,12,12,12,12,12,12,12,12 \\\n'...
        '\t13,13,13,13,13,13,13,13,13,13 \\\n'...
        '\t14,14,14,14,14,14,14,14,14,14 \\\n'...
        '\t15,15,15,15,15,15,15,15,15,15 \\\n'...
        '\t16,16,16,16,16,16,16,16,16,16 \\\n'...
        '\t17,17,17,17,17,17,17,17,17,17 \\\n'...
        '\t43,43,43,43,43,43,43,43,43,43 \\\n'...
        '\t44,44,44,44,44,44,44,44,44,44 \\\n'...
        '\t45,45,45,45,45,45,45,45,45,45 \\\n'...
        '\t46,46,46,46,46,46,46,46,46,46 \\\n'...
        '\t47,47,47,47,47,47,47,47,47,47 \\\n'...
        '\t48,48,48,48,48,48,48,48,48,48 \\\n'...
        '\t49,49,49,49,49,49,49,49,49,49 \\\n'...
        '\t50,50,50,50,50,50,50,50,50,50 \\\n'...
        ')\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'endcue/end_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "end_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.endcue ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sendcue  %s',ttestdir,glmlist_name);
fclose(fid);


%%%% open up Glass ANOVA file 5
% begin script -- Glass
fname = sprintf('%sGlass/glass_fidl_anova.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -tc 26 27 28 29 30 31 32 33  18 19 20 21 22 23 24 25 -scratchdir %s \n\n',glmlist{i},scratchDir.Glass);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t26,26,26,26,26,26,26,26,26,26 \\\n'...
        '\t27,27,27,27,27,27,27,27,27,27 \\\n'...
        '\t28,28,28,28,28,28,28,28,28,28 \\\n'...
        '\t29,29,29,29,29,29,29,29,29,29 \\\n'...
        '\t30,30,30,30,30,30,30,30,30,30 \\\n'...
        '\t31,31,31,31,31,31,31,31,31,31 \\\n'...
        '\t32,32,32,32,32,32,32,32,32,32 \\\n'...
        '\t33,33,33,33,33,33,33,33,33,33 \\\n'...
        '\t18,18,18,18,18,18,18,18,18,18 \\\n'...
        '\t19,19,19,19,19,19,19,19,19,19 \\\n'...
        '\t20,20,20,20,20,20,20,20,20,20 \\\n'...
        '\t21,21,21,21,21,21,21,21,21,21 \\\n'...
        '\t22,22,22,22,22,22,22,22,22,22 \\\n'...
        '\t23,23,23,23,23,23,23,23,23,23 \\\n'...
        '\t24,24,24,24,24,24,24,24,24,24 \\\n'...
        '\t25,25,25,25,25,25,25,25,25,25 \\\n'...
        '\t)\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'Glass/glass_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "glass_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.Glass ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sGlass  %s',ttestdir,glmlist_name);
fclose(fid);


%%%% open up NV ANOVA file 6
% begin script -- NV
fname = sprintf('%sNV/word_fidl_anova.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 51 52 53 54 55 56 57 58  59 60 61 62 63 64 65 66 -scratchdir %s \n\n',glmlist{i},scratchDir.NV);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t51,51,51,51,51,51,51,51,51,51 \\\n'...
        '\t52,52,52,52,52,52,52,52,52,52 \\\n'...
        '\t53,53,53,53,53,53,53,53,53,53 \\\n'...
        '\t54,54,54,54,54,54,54,54,54,54 \\\n'...
        '\t55,55,55,55,55,55,55,55,55,55 \\\n'...
        '\t56,56,56,56,56,56,56,56,56,56 \\\n'...
        '\t57,57,57,57,57,57,57,57,57,57 \\\n'...
        '\t58,58,58,58,58,58,58,58,58,58 \\\n'...
        '\t59,59,59,59,59,59,59,59,59,59 \\\n'...
        '\t60,60,60,60,60,60,60,60,60,60 \\\n'...
        '\t61,61,61,61,61,61,61,61,61,61 \\\n'...
        '\t62,62,62,62,62,62,62,62,62,62 \\\n'...
        '\t63,63,63,63,63,63,63,63,63,63 \\\n'...
        '\t64,64,64,64,64,64,64,64,64,64 \\\n'...
        '\t65,65,65,65,65,65,65,65,65,65 \\\n'...
        '\t66,66,66,66,66,66,66,66,66,66 \\\n'...
        ')\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'NV/word_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "word_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.NV ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sNV  %s',ttestdir,glmlist_name);
fclose(fid);


%%%% open up taskXtime ANOVA file 7
% begin script -- taskXtime
fname = sprintf('%staskXtime/task_fidl_anova.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 18+26 19+27 20+28 21+29 22+30 23+31 24+32 25+33  51+59 52+60 53+61 54+62 55+63 56+64 57+65 58+66 -scratchdir %s \n\n',glmlist{i},scratchDir.tasktime);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t18+26,18+26,18+26,18+26,18+26,18+26,18+26,18+26,18+26,18+26 \\\n'...
        '\t19+27,19+27,19+27,19+27,19+27,19+27,19+27,19+27,19+27,19+27 \\\n'...
        '\t20+28,20+28,20+28,20+28,20+28,20+28,20+28,20+28,20+28,20+28 \\\n'...
        '\t21+29,21+29,21+29,21+29,21+29,21+29,21+29,21+29,21+29,21+29 \\\n'...
        '\t22+30,22+30,22+30,22+30,22+30,22+30,22+30,22+30,22+30,22+30 \\\n'...
        '\t23+31,23+31,23+31,23+31,23+31,23+31,23+31,23+31,23+31,23+31 \\\n'...
        '\t24+32,24+32,24+32,24+32,24+32,24+32,24+32,24+32,24+32,24+32 \\\n'...
        '\t25+33,25+33,25+33,25+33,25+33,25+33,25+33,25+33,25+33,25+33 \\\n'...
        '\t51+59,51+59,51+59,51+59,51+59,51+59,51+59,51+59,51+59,51+59 \\\n'...
        '\t52+60,52+60,52+60,52+60,52+60,52+60,52+60,52+60,52+60,52+60 \\\n'...
        '\t53+61,53+61,53+61,53+61,53+61,53+61,53+61,53+61,53+61,53+61 \\\n'...
        '\t54+62,54+62,54+62,54+62,54+62,54+62,54+62,54+62,54+62,54+62 \\\n'...
        '\t55+63,55+63,55+63,55+63,55+63,55+63,55+63,55+63,55+63,55+63 \\\n'...
        '\t56+64,56+64,56+64,56+64,56+64,56+64,56+64,56+64,56+64,56+64 \\\n'...
        '\t57+65,57+65,57+65,57+65,57+65,57+65,57+65,57+65,57+65,57+65 \\\n'...
        '\t58+66,58+66,58+66,58+66,58+66,58+66,58+66,58+66,58+66,58+66 \\\n'...
        '\t)\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'taskXtime/task_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "task_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.tasktime '  -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%staskXtime  %s',ttestdir,glmlist_name);
fclose(fid);

%%%% open up 4 stimuli ANOVA file 8
% begin script -- ENDCUE
fname = sprintf([ttestdir '4stimXtime/4stim_fidl_anova.csh'],ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 26 27 28 29 30 31 32 33  18 19 20 21 22 23 24 25  51 52 53 54 55 56 57 58  59 60 61 62 63 64 65 66 -scratchdir %s \n\n',glmlist{i},scratchDir.Fourstim);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
outstr = ['set TIME_COURSES = (-tc \\\n'...
        '\t10,10,10,10,10,10,10,10,10,10 \\\n'...
        '\t26,26,26,26,26,26,26,26,26,26 \\\n'...
        '\t27,27,27,27,27,27,27,27,27,27 \\\n'...
        '\t28,28,28,28,28,28,28,28,28,28 \\\n'...
        '\t29,29,29,29,29,29,29,29,29,29 \\\n'...
        '\t30,30,30,30,30,30,30,30,30,30 \\\n'...
        '\t31,31,31,31,31,31,31,31,31,31 \\\n'...
        '\t32,32,32,32,32,32,32,32,32,32 \\\n'...
        '\t33,33,33,33,33,33,33,33,33,33 \\\n'...
        '\t18,18,18,18,18,18,18,18,18,18 \\\n'...
        '\t19,19,19,19,19,19,19,19,19,19 \\\n'...
        '\t20,20,20,20,20,20,20,20,20,20 \\\n'...
        '\t21,21,21,21,21,21,21,21,21,21 \\\n'...
        '\t22,22,22,22,22,22,22,22,22,22 \\\n'...
        '\t23,23,23,23,23,23,23,23,23,23 \\\n'...
        '\t24,24,24,24,24,24,24,24,24,24 \\\n'...
        '\t25,25,25,25,25,25,25,25,25,25 \\\n'...
        '\t51,51,51,51,51,51,51,51,51,51 \\\n'...
        '\t52,52,52,52,52,52,52,52,52,52 \\\n'...
        '\t53,53,53,53,53,53,53,53,53,53 \\\n'...
        '\t54,54,54,54,54,54,54,54,54,54 \\\n'...
        '\t55,55,55,55,55,55,55,55,55,55 \\\n'...
        '\t56,56,56,56,56,56,56,56,56,56 \\\n'...
        '\t57,57,57,57,57,57,57,57,57,57 \\\n'...
        '\t58,58,58,58,58,58,58,58,58,58 \\\n'...
        '\t59,59,59,59,59,59,59,59,59,59 \\\n'...
        '\t60,60,60,60,60,60,60,60,60,60 \\\n'...
        '\t61,61,61,61,61,61,61,61,61,61 \\\n'...
        '\t62,62,62,62,62,62,62,62,62,62 \\\n'...
        '\t63,63,63,63,63,63,63,63,63,63 \\\n'...
        '\t64,64,64,64,64,64,64,64,64,64 \\\n'...
        '\t65,65,65,65,65,65,65,65,65,65 \\\n'...
        '\t66,66,66,66,66,66,66,66,66,66 \\\n'...
        ')\n\n'];
fprintf(fid,outstr);
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir '4stimXtime/4stim_driver.dat" '...
    '-voxel_threshold 0.01 '...
    '-output Z_uncorrected   $GLM_FILES $TIME_COURSES  -Nimage_name "4stim_Nimage.4dfp.img" '...
    '-scratchdir ' scratchDir.Fourstim '  -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 '...
    '-clean_up\n\n'];
fprintf(fid,outstr);
%fprintf(fid,'#%sglm_list.list',ttestdir);
fprintf(fid,'#%s4stimXtime  %s',ttestdir,glmlist_name);

fclose(fid);


end

function write_contrast_script_mem(glmlist_name,procdir,ttestdir,subid,vclist)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% make scratch directories if they do not exist
scratchDir.face = [ttestdir 'face/SCRATCH_face/'];
scratchDir.scene = [ttestdir 'scene/SCRATCH_scene/'];
scratchDir.word = [ttestdir 'word/SCRATCH_word/'];
scratchDir.taskXpres = [ttestdir 'taskXpres/SCRATCH_taskXpres/'];
scratchDir.AllCond = [ttestdir 'AllCond/SCRATCH_AllCond/'];
if ~exist(scratchDir.face)
    mkdir(scratchDir.face);
    mkdir(scratchDir.scene);
    mkdir(scratchDir.word);
    mkdir(scratchDir.taskXpres);
    mkdir(scratchDir.AllCond);
end

% write datfiles
write_datfile_mem(ttestdir,scratchDir,subid,vclist); %need vclist?

%%% Contrasts

% AllCond (vs. Baseline) computation
fname = sprintf('%sAllCond/compute_ttest_single_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir %s -tc  1+2+3+4+5+6+7+8+9+10+11+12+13+14+15+16+17+18+19+20+21+22+23+24+25+26+27+28+29+30+31+32+33+34+35+36+37+38+39+40+41+42+43+44+45+46+47+48+49+50+51+52+53+54+55+56+57+58+59+60+61+62+63+64+65+66+67+68+69+70+71+72+73+74+75+76+77+78+79+80+81+82+83+84+85+86+87+88+89+90+91+92+93+94+95+96+97+98+99+100+101+102+103+104+105+106+107+108+109+110+111+112+113+114+115+116+117+118+119+120+121+122+123+124+125+126+127+128+129+130+131+132+133+134+135+136+137+138+139+140+141+142+143+144\n\n',glmlist{i},scratchDir.AllCond);
end
outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'AllCond/fidl_ttest_single.dat" '...
    '-output Z_uncorrected  '...
    '-Nimage_name "' ttestdir 'AllCond/fidl_ttest_single_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir ' scratchDir.AllCond ' -var_thresh "1e-10" '...
    ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sAllCond/%s',ttestdir,glmlist_name);
fclose(fid);

% ANOVA1: face X pres contrasts
fname = sprintf('%sface/compute_face_anova_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 73 74 75 76 77 78 79 80  81 82 83 84 85 86 87 88  89 90 91 92 93 94 95 96  49 50 51 52 53 54 55 56  57 58 59 60 61 62 63 64  65 66 67 68 69 70 71 72 -scratchdir %s\n\n',glmlist{i},scratchDir.face);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
vals = [73:96 49:72];
fprintf(fid,'set TIME_COURSES = (-tc \\\n');
for v = vals
    fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \\\n',v,v,v,v,v,v,v,v,v,v);
end
fprintf(fid,'\t)\n\n');
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'face/fidl_face.dat" '...
    '-voxel_threshold 0.01 -output Z_uncorrected  $GLM_FILES $TIME_COURSES '...
    '-Nimage_name "faceXpres_Nimage.4dfp.img" -scratchdir ' scratchDir.face ...
    ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%s',glmlist_name);
fclose(fid);

% ANOVA2: scene X pres contrasts
fname = sprintf('%sscene/compute_scene_anova_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 1 2 3 4 5 6 7 8  9 10 11 12 13 14 15 16  17 18 19 20 21 22 23 24  25 26 27 28 29 30 31 32  33 34 35 36 37 38 39 40  41 42 43 44 45 46 47 48 -scratchdir %s\n\n',glmlist{i},scratchDir.scene);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
vals = [1:48];
fprintf(fid,'set TIME_COURSES = (-tc \\\n');
for v = vals
    fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \\\n',v,v,v,v,v,v,v,v,v,v);
end
fprintf(fid,'\t)\n\n');
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'scene/fidl_scene.dat" '...
    '-voxel_threshold 0.01 -output Z_uncorrected  $GLM_FILES $TIME_COURSES '...
    '-Nimage_name "sceneXpres_Nimage.4dfp.img" -scratchdir ' scratchDir.scene ...
    ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%s',glmlist_name);
fclose(fid);

% ANOVA3: word X pres contrasts
fname = sprintf('%sword/compute_word_anova_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 97 98 99 100 101 102 103 104  105 106 107 108 109 110 111 112  113 114 115 116 117 118 119 120  121 122 123 124 125 126 127 128  129 130 131 132 133 134 135 136  137 138 139 140 141 142 143 144 -scratchdir %s\n\n',glmlist{i},scratchDir.word);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
vals = [97:144];
fprintf(fid,'set TIME_COURSES = (-tc \\\n');
for v = vals
    fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d \\\n',v,v,v,v,v,v,v,v,v,v);
end
fprintf(fid,'\t)\n\n');
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'word/fidl_word.dat" '...
    '-voxel_threshold 0.01 -output Z_uncorrected  $GLM_FILES $TIME_COURSES '...
    '-Nimage_name "wordXpres_Nimage.4dfp.img" -scratchdir ' scratchDir.word ...
    ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%s',glmlist_name);
fclose(fid);

% ANOVA4: task X pres contrasts
fname = sprintf('%staskXpres/compute_taskXpres_anova_fidl.csh',ttestdir);
fid = fopen(fname,'w');
outstr = startstring;
fprintf(fid,outstr);
for i = 1:10
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -tc 49+73 50+74 51+75 52+76 53+77 54+78 55+79 56+80  57+81 58+82 59+83 60+84 61+85 62+86 63+87 64+88  65+89 66+90 67+91 68+92 69+93 70+94 71+95 72+96  1+25 2+26 3+27 4+28 5+29 6+30 7+31 8+32  9+33 10+34 11+35 12+36 13+37 14+38 15+39 16+40  17+41 18+42 19+43 20+44 21+45 22+46 23+47 24+48  97+121 98+122 99+123 100+124 101+125 102+126 103+127 104+128  105+129 106+130 107+131 108+132 109+133 110+134 111+135 112+136  113+137 114+138 115+139 116+140 117+141 118+142 119+143 120+144 -scratchdir %s\n\n',glmlist{i},scratchDir.taskXpres);
end
outstr = ['set GLM_FILES = ( -glm_files \\\n'];
fprintf(fid,outstr);
for i = 1:9
    fprintf(fid,'\t %s \\\n',glmlist{i});
end
fprintf(fid,'\t %s)\n\n',glmlist{10});
vals = {'49+73','50+74','51+75','52+76','53+77','54+78','55+79','56+80','57+81','58+82','59+83','60+84','61+85','62+86','63+87','64+88','65+89','66+90','67+91','68+92','69+93','70+94','71+95','72+96','1+25','2+26','3+27','4+28','5+29','6+30','7+31','8+32','9+33','10+34','11+35','12+36','13+37','14+38','15+39','16+40','17+41','18+42','19+43','20+44','21+45','22+46','23+47','24+48','97+121','98+122','99+123','100+124','101+125','102+126','103+127','104+128','105+129','106+130','107+131','108+132','109+133','110+134','111+135','112+136','113+137','114+138','115+139','116+140','117+141','118+142','119+143','120+144'};
fprintf(fid,'set TIME_COURSES = (-tc \\\n');
for v = 1:length(vals)
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s \\\n',vals{v},vals{v},vals{v},vals{v},vals{v},vals{v},vals{v},vals{v},vals{v},vals{v});
end
fprintf(fid,'\t)\n\n');
outstr = ['nice +19 $BIN/fidl_anova4 -driver "' ttestdir 'taskXpres/fidl_taskXpres.dat" '...
    '-voxel_threshold 0.01 -output Z_uncorrected  $GLM_FILES $TIME_COURSES '...
    '-Nimage_name "taskXpres_Nimage.4dfp.img" -scratchdir ' scratchDir.taskXpres ...
    ' -GIGAdesign -glmpersub 1 1 1 1 1 1 1 1 1 1 -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%s',glmlist_name);
fclose(fid);

end

%function write_ttest_script_mixed_OLD(glmlist_name,procdir,ttestdir,subid)

% first load information from glm list
% glm_all = importdata(glmlist_name);
% glmlist = glm_all.textdata(3:12);
% 
% % make a scratch dir if it does not exist 
% scratchDir1 = [procdir 'SCRATCH_ttest1/'];
% scratchDir2 = [procdir 'SCRATCH_ttest2/'];
% if ~exist(scratchDir1)
%     mkdir(scratchDir1);
%     mkdir(scratchDir2);
% end
% 
% % write datfile
% write_datfile_mixed(procdir,scratchDir1,scratchDir2,subid); 
% 
% % open up file
% fname = sprintf([procdir 'compute_ttest_fidl1.csh'],procdir);
% fid = fopen(fname,'w');
% 
% % begin script
% outstr = ['#!/bin/csh\n'...
%     'unlimit\n\n'...
%     'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
%     'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
%     'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
%     'if(`uname` == Linux) then\n'...
%     '\tset dog = `uname -a`\n'...
%     '\tset cat = `expr $#dog - 1`\n'...
%     '\tif($dog[$cat] == x86_64) then\n'...
%     '\t\tset BIN = $BINLINUX64\n'...
%     '\telse\n'...
%     '\t\tset BIN = $BINLINUX\n'...
%     '\tendif\n'...
%     'endif\n\n'...
%     'if($#argv != 1) then\n\n'];
% fprintf(fid,outstr);
% 
% for i = 1:10
%     fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -scratchdir SCRATCH_ttest1/ -tc  20+21+28+29 20+21+28+29+53+54+61+62 3+4 3+4+36+37 36+37 42 53+54+61+62 9 9+42\n\n',glmlist{i});
% end
% 
% 
% outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest1.dat" '...
%     '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_Nimage1.4dfp.img" '...
%     '-glm_list_file "' procdir 'fidl_ttest1.list" -scratchdir SCRATCH_ttest1/ -var_thresh "1e-10" '...
%     '-prepend ' ttestdir ' -glm ' glmlist{i} ' -clean_up\n\n'];
% fprintf(fid,outstr);
% fprintf(fid,'else\n\n');
% fprintf(fid,'nice +19 $BIN/fidl_ttest -driver "%sfidl_ttest1.dat" -clean_up ONLY\n\n',procdir);
% fprintf(fid,'endif\n');
% fprintf(fid,'#%sglm_list.list',procdir);
% 
% fclose(fid);
% 
% % open up file
% fname = sprintf([procdir 'compute_ttest_fidl2.csh'],procdir);
% fid = fopen(fname,'w');
% 
% % begin script
% outstr = ['#!/bin/csh\n'...
%     'unlimit\n\n'...
%     'set BIN = /home/usr/fidl/fidl_code/fidl_2.65/bin\n'...
%     'set BINLINUX = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux\n'...
%     'set BINLINUX64 = /home/usr/fidl/fidl_code/fidl_2.65/bin_linux64\n'...
%     'if(`uname` == Linux) then\n'...
%     '\tset dog = `uname -a`\n'...
%     '\tset cat = `expr $#dog - 1`\n'...
%     '\tif($dog[$cat] == x86_64) then\n'...
%     '\t\tset BIN = $BINLINUX64\n'...
%     '\telse\n'...
%     '\t\tset BIN = $BINLINUX\n'...
%     '\tendif\n'...
%     'endif\n\n'...
%     'if($#argv != 1) then\n\n'];
% fprintf(fid,outstr);
% 
% for i = 1:10
%     fprintf(fid,'nice +19 $BIN/fidl_zstat -glm_file %s -scratchdir SCRATCH_ttest2/ -tc  42 9\n\n',glmlist{i});
% end
% 
% 
% outstr = ['nice +19 $BIN/fidl_ttest -driver "' procdir 'fidl_ttest2.dat" '...
%     '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_Nimage2.4dfp.img" '...
%     '-glm_list_file "' procdir 'fidl_ttest2.list" -scratchdir SCRATCH_ttest2/ -var_thresh "1e-10" '...
%     '-prepend ' ttestdir ' -glm ' glmlist{i} ' -clean_up\n\n'];
% fprintf(fid,outstr);
% fprintf(fid,'else\n\n');
% fprintf(fid,'nice +19 $BIN/fidl_ttest -driver "%sfidl_ttest2.dat" -clean_up ONLY\n\n',procdir);
% fprintf(fid,'endif\n');
% fprintf(fid,'#%sglm_list.list',procdir);
% 
% fclose(fid);

%end


function write_ttest_script_motor(glmlist_name,procdir,ttestdir,subid,vclist)

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);

% make a scratch dir if it does not exist
scratchDir_paired = [ttestdir 'SCRATCH_ttest_paired/'];
scratchDir_single = [ttestdir 'SCRATCH_ttest_single/'];
if ~exist(scratchDir_paired)
    mkdir(scratchDir_paired);
    mkdir(scratchDir_single);
end

% write datfile
write_datfile_motor(ttestdir,scratchDir_paired,scratchDir_single,subid,vclist);

% open up file
fname = [ttestdir 'compute_ttest_paired_fidl.csh'];
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
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir SCRATCH_ttest_paired/ -tc  1 2 2+3 2+3+4+5 2+4 3 3+5 4 4+5 5\n\n',glmlist{i});
end

outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'fidl_ttest_paired.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_paired_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir SCRATCH_ttest_paired/ -var_thresh "1e-10" '...
    '-prepend ' ttestdir ' -clean_up\n\n'];
fprintf(fid,outstr);
fprintf(fid,'#%sglm_list.list',procdir);

fclose(fid);

%%%% begin script -- SINGLE
fname = sprintf([ttestdir 'compute_ttest_single_fidl.csh'],procdir);
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
    fprintf(fid,'nice +19 $BIN/fidl_zstat2 -glm_file %s -scratchdir SCRATCH_ttest_single/ -tc  1 1+2+3+4+5 2 2+3 2+3+4+5 2+4 3 3+5 4 4+5 5\n\n',glmlist{i});
end

outstr = ['nice +19 $BIN/fidl_ttest -driver "' ttestdir 'fidl_ttest_single.dat" '...
    '-output T_uncorrected Z_uncorrected  -Nimage_name "fidl_ttest_single_Nimage.4dfp.img" '...
    '-glm ' glmlist{1} ' -scratchdir SCRATCH_ttest_single/ -var_thresh "1e-10" '...
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

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 0_Tongue+1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_motor_%s_0_Tongue+1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg_1.4dfp.img\n',scratchDir_single,subid,vclist{i});
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

function outstr = startstring
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
end

function write_datfile_mixed(procdir,scratchDir,subid,vclist)

%%%% PAIRED VERSION 1
% open up file
fname = sprintf([procdir 'ttests/fidl_ttest_paired.dat'],procdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'TTEST  := PAIRED_COMPARISON 1_Glass_sustained-7_NV_sustained\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mixed_%s_1_Glass_sustained_1.4dfp.img\n',scratchDir.paired,subid,vclist{i});
end
for i = 1:10
    fprintf(fid,'SECOND := %s%s_mixed_%s_7_NV_sustained_1.4dfp.img\n',scratchDir.paired,subid,vclist{i});
end
fclose(fid);


%%%% SINGLE VERSION 2
% open up file
fname = sprintf([procdir 'ttests/fidl_ttest_single.dat'],procdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_Glass_sustained\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mixed_%s_1_Glass_sustained_1.4dfp.img\n',scratchDir.single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 7_NV_sustained\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mixed_%s_7_NV_sustained_1.4dfp.img\n',scratchDir.single,subid,vclist{i});
end

fprintf(fid,'TTEST  := UNPAIRED_COMPARISON 1_Glass_sustained+7_NV_sustained\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mixed_%s_1_Glass_sustained+7_NV_sustained_1.4dfp.img\n',scratchDir.single,subid,vclist{i});
end
fclose(fid);

%%%% SINGLE VERSION 2b - AllCond vs. baseline
% open up file
fname = sprintf([procdir 'AllCond/fidl_ttest_single.dat'],procdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'TTEST  := UNPAIRED_COMPARISON AllCond\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mixed_%s_toolong1_1+2+3+4+5+6+7+8.4dfp.img\n',scratchDir.AllCond,subid,vclist{i});
end


%%%% STARTCUE ANOVA 3
% open up file
fname = sprintf([procdir 'startcue/start_driver.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'subject	cue	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	glass	%i		%s%s_mixed_%s_0_Glass_startcue_%i.4dfp.img\n',i,j,scratchDir.startcue,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	%i		%s%s_mixed_%s_6_NV_startcue_%i.4dfp.img\n',i,j,scratchDir.startcue,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% ENDCUE ANOVA 4
% open up file
fname = sprintf([procdir 'endcue/end_driver.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'subject	cue	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	glass	%i		%s%s_mixed_%s_2_Glass_endcue_%i.4dfp.img\n',i,j,scratchDir.endcue,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	%i		%s%s_mixed_%s_8_NV_endcue_%i.4dfp.img\n',i,j,scratchDir.endcue,subid,vclist{i},j);
    end
end
fclose(fid);


%%%% GLASS ANOVA 5
% open up file
fname = sprintf([procdir 'Glass/glass_driver.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'subject	stimulus	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	random	%i		%s%s_mixed_%s_4_Glass_random_%i.4dfp.img\n',i,j,scratchDir.Glass,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	coherent	%i		%s%s_mixed_%s_3_Glass_coherent_%i.4dfp.img\n',i,j,scratchDir.Glass,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% NV ANOVA 6
% open up file
fname = sprintf([procdir 'NV/word_driver.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'subject	stimulus	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	noun	%i		%s%s_mixed_%s_9_Noun_%i.4dfp.img\n',i,j,scratchDir.NV,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	verb	%i		%s%s_mixed_%s_10_Verb_%i.4dfp.img\n',i,j,scratchDir.NV,subid,vclist{i},j);
    end
end
fclose(fid);


%%%% taskXtime ANOVA 7
% open up file
fname = sprintf([procdir 'taskXtime/task_driver.dat'],procdir);
fid = fopen(fname,'w');

% write file
fprintf(fid,'subject	task	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	glass	%i		%s%s_mixed_%s_3_Glass_coherent+4_Glass_random_%i.4dfp.img\n',i,j,scratchDir.tasktime,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	%i		%s%s_mixed_%s_9_Noun+10_Verb_%i.4dfp.img\n',i,j,scratchDir.tasktime,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% 4 stimuli ANOVA 8
% open up file
fname = sprintf([procdir '4stimXtime/4stim_driver.dat'],procdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'subject	stimulus	time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	random	%i		%s%s_mixed_%s_4_Glass_random_%i.4dfp.img\n',i,j,scratchDir.Fourstim,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	coherent	%i		%s%s_mixed_%s_3_Glass_coherent_%i.4dfp.img\n',i,j,scratchDir.Fourstim,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	noun	%i		%s%s_mixed_%s_9_Noun_%i.4dfp.img\n',i,j,scratchDir.Fourstim,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	verb	%i		%s%s_mixed_%s_10_Verb_%i.4dfp.img\n',i,j,scratchDir.Fourstim,subid,vclist{i},j);
    end
end
fclose(fid);


end

function write_datfile_mem(ttestdir,scratchDir,subid,vclist)

%%%% FACE
% open up file
fname = sprintf([ttestdir 'face/fidl_face.dat'],ttestdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'subject	face	presentation    time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	male	first	%i	%s%s_mem_%s_9_Male_1st_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	male	second	%i	%s%s_mem_%s_10_Male_2nd_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	male	third	%i	%s%s_mem_%s_11_Male_3rd_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	female	first	%i	%s%s_mem_%s_6_Female_1st_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	female	second	%i	%s%s_mem_%s_7_Female_2nd_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	female	third	%i	%s%s_mem_%s_8_Female_3rd_%i.4dfp.img\n',i,j,scratchDir.face,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% Scene
% open up file
fname = sprintf([ttestdir 'scene/fidl_scene.dat'],ttestdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'subject	scene	presentation    time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	indoor	first	%i	%s%s_mem_%s_0_Indoor_1st_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	indoor	second	%i	%s%s_mem_%s_1_Indoor_2nd_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	indoor	third	%i	%s%s_mem_%s_2_Indoor_3rd_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	outdoor	first	%i	%s%s_mem_%s_3_Outdoor_1st_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	outdoor	second	%i	%s%s_mem_%s_4_Outdoor_2nd_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	outdoor	third	%i	%s%s_mem_%s_5_Outdoor_3rd_%i.4dfp.img\n',i,j,scratchDir.scene,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% Word
% open up file
fname = sprintf([ttestdir 'word/fidl_word.dat'],ttestdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'subject	word	presentation    time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	abstract	first	%i	%s%s_mem_%s_12_Abstract_1st_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	abstract	second	%i	%s%s_mem_%s_13_Abstract_2nd_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	abstract	third	%i	%s%s_mem_%s_14_Abstract_3rd_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	concrete	first	%i	%s%s_mem_%s_15_Concrete_1st_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	concrete	second	%i	%s%s_mem_%s_16_Concrete_2nd_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	concrete	third	%i	%s%s_mem_%s_17_Concrete_3rd_%i.4dfp.img\n',i,j,scratchDir.word,subid,vclist{i},j);
    end
end
fclose(fid);

%%%% SINGLE VERSION 2b - AllCond vs. baseline
% open up file
fname = [ttestdir 'AllCond/fidl_ttest_single.dat'];
fid = fopen(fname,'w');
% write file
fprintf(fid,'TTEST  := UNPAIRED_COMPARISON AllCond\n');
for i = 1:10
    fprintf(fid,'FIRST  := %s%s_mem_%s_toolong1_1+2+3+4+5+6+7+8.4dfp.img\n',scratchDir.AllCond,subid,vclist{i});
end
fclose(fid);

%%%% TaskXPres
% open up file
fname = sprintf([ttestdir 'taskXpres/fidl_taskXpres.dat'],ttestdir);
fid = fopen(fname,'w');
% write file
fprintf(fid,'subject	task	presentation    time	*.4dfp.img\n');
for i = 1:10
    for j = 1:8
        fprintf(fid,'%i	face	first	%i	%s%s_mem_%s_6_Female_1st+9_Male_1st_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	face	second	%i	%s%s_mem_%s_7_Female_2nd+10_Male_2nd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	face	third	%i	%s%s_mem_%s_8_Female_3rd+11_Male_3rd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	scene	first	%i	%s%s_mem_%s_0_Indoor_1st+3_Outdoor_1st_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	scene	second	%i	%s%s_mem_%s_1_Indoor_2nd+4_Outdoor_2nd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	scene	third	%i	%s%s_mem_%s_2_Indoor_3rd+5_Outdoor_3rd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	first	%i	%s%s_mem_%s_12_Abstract_1st+15_Concrete_1st_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	second	%i	%s%s_mem_%s_13_Abstract_2nd+16_Concrete_2nd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
    for j = 1:8
        fprintf(fid,'%i	word	third	%i	%s%s_mem_%s_14_Abstract_3rd+17_Concrete_3rd_%i.4dfp.img\n',i,j,scratchDir.taskXpres,subid,vclist{i},j);
    end
end
fclose(fid);



end
