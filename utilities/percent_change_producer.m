function percent_change_producer()

sub_nums = [2]; %1:10;
tasks = {'motor','mem','mixed'};
condition_sets = {{'AllCond'},{'AllCond','AllFace','AllScene','AllWord'},{'AllCond','AllGlass','AllSemantic'}};

topDir = '/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/';

for s = 1:length(sub_nums)
    
    sub_id = sprintf('MSC%02d',sub_nums(s));
    subDir = [topDir sub_id '/'];
    if ~exist(subDir)
        mkdir(subDir);
    end
    dataDir = ['/data/nil-bluearc/GMT/Laumann/MSC/' sub_id '/'];
    disp(['Sub: ' sub_id]);
    
    for t = 1:length(tasks)
        
        % task directory
        taskDir = [topDir sub_id '/' tasks{t} '_surf/'];
        
        conds = condition_sets{t};
        
        % output directory
        outDir = [taskDir 'percent_change/'];
        mkdir(outDir);
        
        % write script
        glm_fname = [taskDir 'glm_list.list'];
        
        switch tasks{t}
            case 'mem'
                if strcmp(sub_id,'MSC10')
                    disp('MSC10 needs to be done by hand for mem');
                else
                    write_perc_ch_script_mem(glm_fname,outDir,sub_id,tasks{t},conds);
                end
            case 'motor'
                write_perc_ch_script_motor(glm_fname,outDir,sub_id,tasks{t},conds);
            case 'mixed'
                write_perc_ch_script_mixed(glm_fname,outDir,sub_id,tasks{t},conds);
        end
    end
    
    
end


end

function write_perc_ch_script_motor(glmlist_name,outDir,sub,task,conds)
% hard coded for just AllCond condition currently

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);
for i = 1:length(glmlist)
    tmp_ind = strfind(glmlist{i},'vc');
    tmp_ind2 = strfind(glmlist{i},'.glm');
    vcidlist{i} = glmlist{i}(tmp_ind:tmp_ind2-1);
end

% open up file
fname = [outDir 'fidl_avg_zstat_' task '.csh'];
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
    fprintf(fid,'set GLM_FILES = (-glm_files %s)\n',glmlist{i});
    fprintf(fid,'set TIME_COURSES = (-tc 1+2+3+4+5)\n');
    fprintf(fid,'nice +19 $BIN/fidl_avg_zstat2 $GLM_FILES $TIME_COURSES -group_name "" -glm_list_file %s -tc_names "%s_%s_%s_%s" -frames 1 -glmpersub 1\n\n',glmlist_name,sub,task,vcidlist{i},conds{1});
end

fclose(fid);

currDir = pwd;
cd(outDir);
system(['csh ' fname]);
cd(currDir);

end

function write_perc_ch_script_mem(glmlist_name,outDir,sub,task,conds)
%%% looking at tp 3+4 for now as the estimate
% hard coded to do AllCond, AllFace, AllScene, AllWord currently

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);
for i = 1:length(glmlist)
    tmp_ind = strfind(glmlist{i},'vc');
    tmp_ind2 = strfind(glmlist{i},'.glm');
    vcidlist{i} = glmlist{i}(tmp_ind:tmp_ind2-1);
end
    
% open up file
fname = [outDir 'fidl_avg_zstat_' task '.csh'];
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
    fprintf(fid,'set GLM_FILES = (-glm_files %s)\n',glmlist{i});
    fprintf(fid,'set TIME_COURSES = (-tc \\\n');
    fprintf(fid,'\t3+4+11+12+19+20+27+28+35+36+43+44+51+52+59+60+67+68+75+76+83+84+91+92+99+100+107+108+115+116+123+124+131+132+139+140 \\\n');
    fprintf(fid,'\t3+4+11+12+19+20+27+28+35+36+43+44 \\\n');
    fprintf(fid,'\t51+52+59+60+67+68+75+76+83+84+91+92 \\\n');
    fprintf(fid,'\t99+100+107+108+115+116+123+124+131+132+139+140)\n');
    fprintf(fid,'nice +19 $BIN/fidl_avg_zstat2 $GLM_FILES $TIME_COURSES -group_name "" -glm_list_file %s -tc_names "%s_%s_%s_%s" "%s_%s_%s_%s" "%s_%s_%s_%s" "%s_%s_%s_%s" -frames 1 1 1 1 -glmpersub 1\n\n',...
        glmlist_name,sub,task,vcidlist{i},conds{1},...
        sub,task,vcidlist{i},conds{2},...
        sub,task,vcidlist{i},conds{3},...
        sub,task,vcidlist{i},conds{4});
end

fclose(fid);

currDir = pwd;
cd(outDir);
system(['csh ' fname]);
cd(currDir);
end

function write_perc_ch_script_mixed(glmlist_name,outDir,sub,task,conds)
%%% looking at tp 3+4 for now as the estimate
% hard coded to do AllCond, AllGlass, and AllSemantic currently

% first load information from glm list
glm_all = importdata(glmlist_name);
glmlist = glm_all.textdata(3:12);
for i = 1:length(glmlist)
    tmp_ind = strfind(glmlist{i},'vc');
    tmp_ind2 = strfind(glmlist{i},'.glm');
    vcidlist{i} = glmlist{i}(tmp_ind:tmp_ind2-1);
end

% open up file
fname = [outDir 'fidl_avg_zstat_' task '.csh'];
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
    fprintf(fid,'set GLM_FILES = (-glm_files %s)\n',glmlist{i});
    fprintf(fid,'set TIME_COURSES = (-tc \\\n');
    fprintf(fid,'\t3+4+11+12+13+20+21+28+29+36+37+44+45+46+53+54+61+62 \\\n');
    fprintf(fid,'\t3+4+11+12+13+20+21+28+29 \\\n');
    fprintf(fid,'\t36+37+44+45+46+53+54+61+62)\n');
    fprintf(fid,'nice +19 $BIN/fidl_avg_zstat2 $GLM_FILES $TIME_COURSES -group_name "" -glm_list_file %s -tc_names "%s_%s_%s_%s" "%s_%s_%s_%s" "%s_%s_%s_%s" -frames 1 1 1 -glmpersub 1\n\n',...
        glmlist_name,sub,task,vcidlist{i},conds{1},...
        sub,task,vcidlist{i},conds{2},...
        sub,task,vcidlist{i},conds{3});
end

fclose(fid);

currDir = pwd;
cd(outDir);
system(['csh ' fname]);
cd(currDir);

end