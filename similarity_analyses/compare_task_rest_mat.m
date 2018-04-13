function compare_task_rest_mat(task,matchType,roiType)
% function compare_task_rest_mat(task,matchType,roiType)
%%% Quick script to compare task and rest matrices
% saves out correlation matrices needed for later comparisons
% including correlation matrices concatenated across all sessions, specific
% sessions, and for split halves
%
% Inputs:
%   task = string with task name ('motor','mixed', or 'mem')
%   matchType = string with 'orig' or 'match' (full dataset or matched frames)
%   roiType = '333', '333_preGLM' or 'Ind' (run correlations on 333 group parcels or individual parcels)
%
% C. Gratton

% subjects to run analyses on
subjects = [1:10];

% some constants
TR = 2.2;
start_scan_cut = round(30/TR); % number of frames to cut from the start of each scan to account for the start scan effect
splithalf1_sessions = [1:5]; %sessions going into each split half
splithalf2_sessions = [6:10];

% directory structure
topDir = '/data/nil-bluearc/GMT/';
switch roiType
    case '333'
        FCDir = [topDir 'Caterina/TaskFC/FC_Parcels/'];
    case 'Ind'
        FCDir = [topDir 'Caterina/TaskFC/FC_Parcels_Ind/'];
    case '333_preGLM'
        FCDir = [topDir 'Caterina/TaskFC/FC_Parcels_preGLM/'];
end

switch task
    case 'motor'
        cond_types_all = {'rest','AllMotor','Tongue','L_Hand','R_Hand','L_Leg','R_Leg'};
        cond_types = {'AllMotor','Tongue','L_Hand','R_Hand','L_Leg','R_Leg'};
    case 'mixed'
        cond_types_all = {'rest','AllGlass','AllSemantic'};
        cond_types = {'AllGlass','AllSemantic'};
    case 'mem'
        cond_types_all = {'rest','AllMem','Scene','Face','Word','pres1','pres2','pres3'};
        cond_types = {'AllMem','Scene','Face','Word','pres1','pres2','pres3'};
end

% reduce subject set and condition set if match type
if strcmp(matchType,'match')
    disp('good match subs, >24.5 min. of data:'); 
    %%%
    match_fr = 675; %700; %decreased to let MSC07 stay in
    % other bined amts will take match_fr/nbins*0.9
    %match_fr_half = round(match_fr/2)-50; % give a little more leeway for halves to be different size
    subjects = [1:7,10]
    if strcmp(task,'motor')
        cond_types = {'AllMotor'};
        cond_types_all = {'rest','AllMotor'};        
    elseif strcmp(task,'mem')
        cond_types = {'AllMem'}'%,'Scene','Face','Word'};
        cond_types_all = {'rest','AllMem'};%,'Scene','Face','Word'};
        % exclude 3 and 9 if I do scene, face, word separately
    end

end

% ROI information
atlas_params = atlas_parameters('Parcels','/data/cn5/caterina/Atlases/Evan_parcellation/');

for s = 1:length(subjects)
    
    subject = sprintf('MSC%02d',subjects(s));
    sub_list{s} = subject;
    disp(['subject: ' subject]);
    
    taskFCDir = [FCDir task '/'];
    restFCDir = [FCDir 'rest/'];

    % load QC file for this subject to compute start scan mask - TASK
    load(['/data/nil-bluearc/GMT/Caterina/TaskFC/FCProc_' subject '_' task '_pass2/QC.mat']);
    start_scan_mask = construct_start_scan_mask(QC,start_scan_cut);
    
    % load session numbering information
    session_numbers = get_session_numbers(subject,QC,task);
    %session_numbers_all{s} = session_numbers; % store for checking later

    % delete QC information file
    clear QC;

    % load task indices
    task_conds = load([topDir 'Caterina/TaskFC/FCProc_' subject '_' task '_pass2/condindices.mat']);
    %cond_types = task_conds.cond_types;
    
    %load task data
    %load([taskFCDir subject '_parcel_corrmat.mat']);
    task_ts = load([taskFCDir subject '_parcel_timecourse.mat']);    
    for c = 1:length(cond_types)  
        disp(['cond: ' cond_types{c}]);
        
        %concatenate ts across sessions
        ts = []; ts_half1 = []; ts_half2 = []; ts_bysess = [];
        FD_half1 = []; FD_half2 = []; %FrameNum_half1 = []; FrameNum_half2 = [];
        count = 1; goodsess_counter = []; session_numbers_final = [];
        for i = 1:length(task_ts.parcel_time) 
            if size(task_ts.parcel_time{i},1)>1
                temp_tmask = task_conds.TIndFin(i).(cond_types{c}).*start_scan_mask{count}; %has combo of task and tmask already in it, add in start scan mask
                ts = [ts; task_ts.parcel_time{i}(logical(temp_tmask),:)];
                %task_ts_bysess.(cond_types{c}){s,i} = task_ts.parcel_time{i}(logical(temp_tmask),:);
                if sum(temp_tmask) > 1 % make sure there is SOME good data in here
                    ts_bysess{i} = task_ts.parcel_time{i}(logical(temp_tmask),:);
                    session_numbers_final = [session_numbers_final session_numbers(count)]; %this is a label for the sess - ie to deal with MSC02
                    goodsess_counter = [goodsess_counter i]; %this is a 1-10 index
                else
                    ts_bysess{i} = []; 
                end
                FD_bysess(i) = task_conds.FDtot(i).(cond_types{c}); % FD of full run, not masked
                %FrameNum_bysess{i} = sum(temp_tmask);
                count = count+1; % counter of sessions analyzed - not necessarily all good post temp_tmask
                                
                % also CREATE a split half version
                %if find([1 4 5 8 9] == i) %approx matched in order effects = ORIGINAL VERSION FOR NEURON SUBMISSION
                if sum(i==splithalf1_sessions)>0 % NEW VERSION, coded at the top
                    ts_half1 = [ts_half1; task_ts.parcel_time{i}(logical(temp_tmask),:)];
                    FD_half1 = [FD_half1; task_conds.FDtot(i).(cond_types{c})];
                elseif sum(i == splithalf2_sessions)>0
                    ts_half2 = [ts_half2; task_ts.parcel_time{i}(logical(temp_tmask),:)];
                    FD_half2 = [FD_half1; task_conds.FDtot(i).(cond_types{c})];
                else
                    error('something is wrong with your session numbering');
                end
            else
                ts_bysess{i} = [];
                FD_bysess(i) = task_conds.FDtot(i).(cond_types{c});
            end
        end
        
        % shorten TS if in match case, try to evently sample by session
        switch matchType
            case 'match'
                disp('match all');
                task_ts_all.(cond_types{c}){s} = match_data(ts_bysess,match_fr,1); %ts(1:match_fr,:);
                disp('match halves');
                %task_ts_half1.(cond_types{c}){s} = match_data(ts_bysess([1 4 5 8 9]),match_fr,2); %ts_half1(1:match_fr_half,:); -> OLD VERSION               
                %task_ts_half2.(cond_types{c}){s} = match_data(ts_bysess([2 3 6 7 10]),match_fr,2); %ts_half2(1:match_fr_half,:); -> OLD VERSION
                task_ts_half1.(cond_types{c}){s} = match_data(ts_bysess(splithalf1_sessions),match_fr,2);          
                task_ts_half2.(cond_types{c}){s} = match_data(ts_bysess(splithalf2_sessions),match_fr,2);
                disp('match sessions');
                [task_ts_bysess.(cond_types{c})(s,:) session_numbers_all.(cond_types{c}){s} goodsessions.(cond_types{c}){s}] = match_data_bysess(ts_bysess,match_fr,10,session_numbers_final,goodsess_counter);
                %if length(goodsessions.(cond_types{c}){s}) ~= length(goodsess_counter)
                %    disp('Check whats up');
                %end                
            case 'orig'
                task_ts_all.(cond_types{c}){s} = ts;
                task_ts_half1.(cond_types{c}){s} = ts_half1;
                task_ts_half2.(cond_types{c}){s} = ts_half2;
                task_ts_bysess.(cond_types{c})(s,:) = ts_bysess;
                session_numbers_all.(cond_types{c}){s} = session_numbers_final; %[1:10];
                goodsessions.(cond_types{c}){s} = goodsess_counter; %[1:10];
        end
        FD_mean_all.(cond_types{c}){s} = mean(FD_bysess(goodsess_counter)); % these don't change with matching, since we care about pre-selection FD
        FD_mean_half1.(cond_types{c}){s} = mean(FD_half1);
        FD_mean_half2.(cond_types{c}){s} = mean(FD_half2);
        FD_mean_bysess.(cond_types{c})(s,:) = FD_bysess;
        
    end   
    
    % calculate the correlations
    for c = 1:length(cond_types)
        avg_task_corrmat.(cond_types{c})(:,:,s) = atanh(paircorr_mod(task_ts_all.(cond_types{c}){s}));
        avg_task_corrmat_halves.(cond_types{c})(:,:,s,1) = atanh(paircorr_mod(task_ts_half1.(cond_types{c}){s}));
        avg_task_corrmat_halves.(cond_types{c})(:,:,s,2) = atanh(paircorr_mod(task_ts_half2.(cond_types{c}){s}));
        for i = 1:10 %length(task_ts.parcel_time) - - ASSUME THERE are 10 sessions
            if find(goodsessions.(cond_types{c}){s}==i)
                sess_task_corrmat.(cond_types{c})(:,:,s,i) = atanh(paircorr_mod(task_ts_bysess.(cond_types{c}){s,i}));
            else %make this all nans
                sess_task_corrmat.(cond_types{c})(:,:,s,i) = ones(size(task_ts_all.(cond_types{c}){s},2))*nan;
            end
        end
    end
    
    % load QC file for this subject to compute start scan mask - REST
    load(['/data/nil-bluearc/GMT/Laumann/MSC/' subject '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/QC.mat']);
    start_scan_mask = construct_start_scan_mask(QC,start_scan_cut);
    FD_bysess = get_FD_vals(QC);
    clear QC;    
    
    % load rest data
    rest_ts = load([restFCDir subject '_parcel_timecourse.mat']);
    ts = []; ts_half1 = []; ts_half2 = [];
    FD_half1 = []; FD_half2 = [];
    for i = 1:length(rest_ts.parcel_time) % concatanate ts across sessions - ASSUME THERE are 10 sessions for REST - no sessions ever lost
        ts = [ts; rest_ts.parcel_time{i}(logical(rest_ts.tmask_all{i}),:)];
        ts_bysess{i} = rest_ts.parcel_time{i}(logical(rest_ts.tmask_all{i}),:);
        
         % also do a split half version
         %if find([1 4 5 8 9] == i) %approx matched in order effects = ORIGINAL VERSION FOR NEURON SUBMISSION
         if sum(i==splithalf1_sessions)>0 % first half vs. second half % NEW VERSION
             ts_half1 = [ts_half1; rest_ts.parcel_time{i}(logical(rest_ts.tmask_all{i}),:)];
             FD_half1 = [FD_half1; FD_bysess(i)];
         elseif sum(i==splithalf2_sessions)>0
             ts_half2 = [ts_half2; rest_ts.parcel_time{i}(logical(rest_ts.tmask_all{i}),:)];
             FD_half2 = [FD_half2; FD_bysess(i)];
         else
             error('something went wrong with session numbering');
         end        
    end    
    
    switch matchType
            case 'match'
                rest_ts_all{s} = match_data(ts_bysess,match_fr,1); %ts(1:match_fr,:);
                %rest_ts_half1{s} = match_data(ts_bysess([1 4 5 8 9]),match_fr,2); %ts_half1(1:match_fr_half,:); = OLD VERSION               
                %rest_ts_half2{s} = match_data(ts_bysess([2 3 6 7 10]),match_fr,2); %ts_half2(1:match_fr_half,:); = OLD VERSION
                rest_ts_half1{s} = match_data(ts_bysess(splithalf1_sessions),match_fr,2); %ts_half1(1:match_fr_half,:);                
                rest_ts_half2{s} = match_data(ts_bysess(splithalf2_sessions),match_fr,2); %ts_half2(1:match_fr_half,:);
                [rest_ts_bysess(s,:) rest_session_numbers_all{s} rest_goodsessions{s}] = match_data_bysess(ts_bysess,match_fr,10,[1:10],[1:10]);
            case 'orig'
                rest_ts_all{s} = ts;
                rest_ts_half1{s} = ts_half1;
                rest_ts_half2{s} = ts_half2;
                rest_ts_bysess(s,:) = ts_bysess;
                rest_goodsessions{s} = [1:10];
                rest_session_numbers_all{s} = [1:10];
    end
    rest_FD_mean_all{s} = mean(FD_bysess); % these don't change with matching, since we care about pre-selection FD
    rest_FD_mean_half1{s} = mean(FD_half1);
    rest_FD_mean_half2{s} = mean(FD_half2);
    rest_FD_mean_bysess(s,:) = FD_bysess;

    
    
    avg_rest_corrmat(:,:,s) = atanh(paircorr_mod(rest_ts_all{s})); %calculate the correlations
    avg_rest_corrmat_halves(:,:,s,1) = atanh(paircorr_mod(rest_ts_half1{s}));
    avg_rest_corrmat_halves(:,:,s,2) = atanh(paircorr_mod(rest_ts_half2{s}));
    for i = 1:10 %length(rest_ts.parcel_time)- ASSUME THERE are 10 sessions
        sess_rest_corrmat(:,:,s,i) = atanh(paircorr_mod(rest_ts_bysess{s,i}));
    end

    
    % plot the original datasets and save    
    %parcel_correlmat_figmaker_cg(squeeze(avg_task_corrmat(s,:,:)),['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    %save_fig(gcf,[taskFCDir '/' subject '_' task '_parcel_corrmat.png']);
    %parcel_correlmat_figmaker_cg(squeeze(avg_rest_corrmat(s,:,:)),['/data/cn5/caterina/TaskConn_Methods/all_data/ParcelCommunities.txt'],[-0.4 1]);
    %save_fig(gcf,[restFCDir_out '/' subject '_rest_parcel_corrmat.png']);
    
    % plot and save the difference
    for c = 1:length(cond_types)
        diffmat = squeeze(avg_task_corrmat.(cond_types{c})(:,:,s)-avg_rest_corrmat(:,:,s));
        figure_corrmat_network_generic(diffmat,atlas_params,[-0.4 0.4]);
        save_fig(gcf,[taskFCDir '/' subject '_' task '_cond' cond_types{c} 'minusrest_parcel_corrmat_' matchType '.png']);
    end
    
    close('all');
    
end

% save datasets
save([restFCDir '/allsubs_rest_corrmats_' matchType '.mat'],'avg_rest_corrmat','rest_goodsessions','rest_session_numbers_all','rest_FD_mean_all');
save([taskFCDir '/allsubs_' task '_corrmats_' matchType '.mat'],'avg_task_corrmat','goodsessions','session_numbers_all','FD_mean_all');

save([restFCDir '/allsubs_rest_corrmats_halves_' matchType '.mat'],'avg_rest_corrmat_halves','rest_ts_half1','rest_ts_half2','rest_goodsessions','rest_session_numbers_all','rest_FD_mean_half1','rest_FD_mean_half2');
save([taskFCDir '/allsubs_' task '_corrmats_halves_' matchType  '.mat'],'avg_task_corrmat_halves','task_ts_half1','task_ts_half2','goodsessions','session_numbers_all','FD_mean_half1','FD_mean_half2');

save([restFCDir '/allsubs_rest_corrmats_bysess_' matchType '.mat'],'sess_rest_corrmat','rest_ts_bysess','rest_goodsessions','rest_session_numbers_all','rest_FD_mean_bysess');
save([taskFCDir '/allsubs_' task '_corrmats_bysess_' matchType '.mat'],'sess_task_corrmat','task_ts_bysess','goodsessions','session_numbers_all','FD_mean_bysess');

    
% organize_data
big_data_mat = avg_rest_corrmat;
big_data_cond_ind = ones(1,size(big_data_mat,3));
big_data_sub_ind = [1:size(big_data_mat,3)];
for c = 1:length(cond_types)
    nsub = size(avg_task_corrmat.(cond_types{c}),3);
    big_data_mat(:,:,end+1:end+nsub) = avg_task_corrmat.(cond_types{c});
    big_data_cond_ind(end+1:end+nsub) = ones(1,nsub).*(c+1);
    big_data_sub_ind(end+1:end+nsub) = [1:nsub];
end
    
% MDS plot
figure('Position',[1 1 1400 600]);
subplot(1,2,1)
cmdscale_mat(big_data_mat,big_data_sub_ind,'euclidean',sub_list);
title('By Subject')
grid on;
subplot(1,2,2)
cmdscale_mat(big_data_mat,big_data_cond_ind,'euclidean',cond_types_all);
title('By Condition')
grid on;
save_fig(gcf,[taskFCDir '/MDS_allconds_' matchType '.pdf']);


% Run OODA on all group comparisons
%run_OODA_weighted_groups(big_data_mat,big_data_cond_ind,['allconds_' task],cond_types_all,...
%    ['/data/nil-bluearc/GMT/Caterina/TaskFC/FC_Parcels/OODA/' task '_allconds_' matchType '/']);

% and do again, but compared with rest
big_data_diff_mat = [];
big_data_diff_cond_ind = [];
big_data_diff_sub_ind = [];
for c = 1:length(cond_types)
    nsub = size(avg_task_corrmat.(cond_types{c}),3);
    big_data_diff_mat(:,:,end+1:end+nsub) = avg_task_corrmat.(cond_types{c}) - avg_rest_corrmat;
    big_data_diff_cond_ind(end+1:end+nsub) = ones(1,nsub).*c;
    big_data_diff_sub_ind(end+1:end+nsub) = [1:nsub];
end
figure('Position',[1 1 1400 600]);
subplot(1,2,1)
cmdscale_mat(big_data_diff_mat,big_data_diff_sub_ind,'euclidean',sub_list);

title('By Subject, diffmats')
grid on;
subplot(1,2,2)
cmdscale_mat(big_data_diff_mat,big_data_diff_cond_ind,'euclidean',cond_types);
grid on;
title('By Condition, diffmats')
save_fig(gcf,[taskFCDir '/MDS_allconds_diffmats_' matchType '.pdf']);

% Run OODA on all group comparisons
%run_OODA_weighted_groups(big_data_diff_mat,big_data_diff_cond_ind,['diffs_' task],cond_types,...
%    ['/data/nil-bluearc/GMT/Caterina/TaskFC/FC_Parcels/OODA/' task '_diffs/']);

end

function ts_concat = match_data(timeseries,totalamt,nbins)

nsess = 10; % always 10

% calculate how much you want total for this binning
bin_amt = floor(totalamt/nbins);
if nbins > 1
    bin_amt = floor(bin_amt * 0.8); % only take 80% to give more leeway in being able to grab full amt of data per sess
end

npersess = floor(bin_amt/(nsess/nbins)); %now grab approximately equal amts from all of the sessions
ts_concat = [];
for s = 1:length(timeseries)
    if size(timeseries{s},1) < npersess
        ts_concat = [ts_concat; timeseries{s}]; % take everything and you'll make up for it in the end
        disp(['too few TRs in this session: ' num2str(s)]);
    else
        ts_concat = [ts_concat; timeseries{s}(1:npersess,:)];
    end
end

%and grab whatever you need (from the first session with enough) to get to full amt
leftover = bin_amt - size(ts_concat,1);
needmore = 1;
while leftover > 0
    if size(timeseries{needmore},1) >= (leftover+npersess);
        ts_concat = [ts_concat; timeseries{needmore}(npersess+1:npersess+leftover,:)];
    else
        ts_concat = [ts_concat; timeseries{needmore}(npersess+1:end,:)];
        needmore = needmore + 1;
    end
    leftover = bin_amt - size(ts_concat,1);
end
end

function [ts_concat session_numbers_final goodsessions_final] = match_data_bysess(timeseries,totalamt,nbins,session_numbers,goodsessions_init)

nsess = length(session_numbers); %10; % always 10

% calculate how much you want total for this binning
bin_amt = floor(totalamt/nbins);
if nbins > 1
    bin_amt = floor(bin_amt * 0.75); % only take 75% to give more leeway in being able to grab full amt of data per sess
end

npersess = floor(bin_amt/(nsess/nbins)); %now grab approximately equal amts from all of the sessions
session_numbers_final = []; goodsessions_final = [];
count = 1;
for s = 1:length(timeseries)
    %deal with cases where there was insufficient data already to analyze
    if size(timeseries{s},1) > 2
        if size(timeseries{s},1) < npersess
            ts_concat{s} = []; %ones(npersess,size(timeseries{s},2)).*nan; % make this a bad session if it doesn't reach sufficient levels
            disp(['too few TRs in this session: ' num2str(s)]);
        else
            ts_concat{s} = timeseries{s}(1:npersess,:);
            session_numbers_final = [session_numbers_final session_numbers(count)]; %actual name
            goodsessions_final = [goodsessions_final s]; %index 1-10
        end
        count = count + 1;
    %else
    %    disp(['Let me know']);
    end
end

end

function [session_numbers ] = get_session_numbers(subject,QC,task)

%this is unnecessary in many cases, but helps when the datalist was cut to
%ensure we're matching across sessions properly
topDir = '/data/nil-bluearc/GMT/Caterina/TaskFC/';
datalist = [topDir subject '_' task '_DATALIST.txt'];
[datapath vcid fcparams TR frames] = textread(datalist,'%s%s%s%f%d');

session_numbers = []; %goodsessions = [];
for i = 1:length(QC)
    ind = find(strcmp(QC(i).vcnum,vcid));
    session_numbers = [session_numbers ind];
end

%%% And correct for weird exceptions from MSC02 (but only for session numbers, not ind):
if strcmp(subject,'MSC02')
    if strcmp(task,'motor')
        session_numbers(session_numbers == 4) = 11; % switch 4 to novel 11 session
    elseif strcmp(task,'mixed') %weird split across two sessions
        session_numbers(session_numbers == 3) = 12;
    elseif strcmp(task,'mem') %weird split across two session numbers
        session_numbers(session_numbers == 8) = 13;
    end
end

end

function FD = get_FD_vals(QC)

for i = 1:length(QC)
    FD(i) = QC(i).FDbar;
end

end
