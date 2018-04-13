function create_task_activation_mats(task)%,roiType)
% function create_task_activation_mats(task)%,roiType)
% Quick script to compare task and rest matrices
% task = 'motor' or 'mixed' or 'mem' 
%
% CG

% subjects
subjects = [1:10]; %[1:10];

% directory information
topDir = '/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/';
allsub_outDir = ['/data/nil-bluearc/GMT/Caterina/TaskFC/FC_Parcels/' task '/'];

% condition information
switch task
    case 'motor'
        cond_types = {'AllCond'}; %{'0_Tongue+1_L_Hand+2_R_Hand+3_L_Leg+4_R_Leg'}; %{'AllCond'};
    case 'mixed'
        cond_types = {'AllCond','AllGlass','AllSemantic'};
    case 'mem'
        cond_types = {'AllCond','AllFace','AllScene','AllWord'};
end

% split half condition information
sh1 = [1:5];
sh2 = [6:10];
%sh1 = [1 4 6 8 9];
%sh2 = [2 3 5 7 10];

for s = 1:length(subjects)
    
    subject = sprintf('MSC%02d',subjects(s));
    sub_list{s} = subject;
    
    % get location of data and order of VCIDs - note that this varies
    % slightly by task in a couple of cases
    [vclist session_numbers] = get_vclist(subject,task);
    
    % where percent signal change files are saved
    taskActDir = [topDir subject '/' task '_surf/percent_change/' ];
    
    %load task data per parcel
    watershed_LR = '/data/cn4/laumannt/Parcellation/Parcels_LR.dtseries.nii';
    startname = [subject '_' task];
    for c = 1:length(cond_types)
        [parcel_act.(cond_types{c}), parcel_act_concat.(cond_types{c})] = extract_activations_parcels(watershed_LR,vclist,taskActDir,startname,cond_types{c});
                
        % also do a split half version
        parcel_act_concat_halves.(cond_types{c})(:,1) = nanmean(parcel_act_concat.(cond_types{c})(sh1,:),1);
        parcel_act_concat_halves.(cond_types{c})(:,2) = nanmean(parcel_act_concat.(cond_types{c})(sh2,:),1);

        % make a master version across subjects
        allsubs_parcel_act.(cond_types{c}){s} = parcel_act.(cond_types{c});
        allsubs_parcel_act_concat.(cond_types{c})(:,s,:) = parcel_act_concat.(cond_types{c})';
        allsubs_parcel_act_concat_halves.(cond_types{c})(:,s,:) = parcel_act_concat_halves.(cond_types{c});

    end
    
    % save datasets
    save([taskActDir subject '_percent_change_groupParcels.mat'],'parcel_act_concat','parcel_act','session_numbers');
    save([taskActDir subject '_percent_change_groupParcels_halves.mat'],'parcel_act_concat_halves','session_numbers');
    
        
end

save([allsub_outDir 'percent_change_allsubs.mat'],'allsubs_parcel_act','allsubs_parcel_act_concat','allsubs_parcel_act_concat_halves');

% organize_data
% big_data_mat = [];
% big_data_cond_ind = [];
% big_data_sub_ind = [];
% for c = 1:length(cond_types)
%     nsub = size(allsubs_parcel_act_concat.(cond_types{c}),2);
%     big_data_mat(:,end+1:end+nsub) = mean(allsubs_parcel_act_concat.(cond_types{c}),3);
%     big_data_cond_ind(end+1:end+nsub) = ones(1,nsub).*c;
%     big_data_sub_ind(end+1:end+nsub) = [1:nsub];
% end
% 
% % MDS plot
% figure('Position',[1 1 1400 600]);
% subplot(1,2,1)
% cmdscale_lin(big_data_mat,big_data_sub_ind,'euclidean',sub_list);
% title('By Subject')
% grid on;
% subplot(1,2,2)
% cmdscale_lin(big_data_mat,big_data_cond_ind,'euclidean',cond_types);
% title('By Condition')
% grid on;
% save_fig(gcf,[taskFCDir '/MDS_allconds_' matchType '.pdf']);


end

function [vclist session_numbers] = get_vclist(sub_id,task)

session_numbers = [1:10]; % generally - minus a couple of weirdnesses noted below

if strcmp(sub_id,'MSC02')
    if strcmp(task,'motor')
        func_file = ['/data/nil-bluearc/GMT/MSC_task/MSC_glms_Becky/' sub_id '/' sub_id '_motor_functional_vclist.txt'];
        session_numbers(session_numbers == 4) = 11;
    else
        func_file = ['/data/nil-bluearc/GMT/MSC_task/MSC_glms_Becky/' sub_id '/' sub_id '_mixed_memory_functional_vclist.txt'];
        if strcmp(task,'mixed')
              session_numbers(session_numbers == 3) = 12;
        elseif strcmp(task,'mem')
              session_numbers(session_numbers == 8) = 13;
        end

    end
else
    func_file = ['/data/nil-bluearc/GMT/MSC_task/MSC_glms_Becky/' sub_id '/' sub_id '_functional_vclist.txt'];
end

[vclist_init] = textread(func_file,'%s');

for v = 1:length(vclist_init)
    vclist{v} = ['vc' vclist_init{v}];
end

end
