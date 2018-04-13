function [allvals] = corr_relationship_bynod(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir_orig,title_start)
% function [allvals] = corr_relationship_bynod(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir_orig,title_start)
% function to calculate the similarity by row (node) and show projected onto each parcel
% Inputs:
%   big_data_mat: a node x node x type matrix, where type = a particular
%      subject-session-task grouping for a dataset (made in similarity_figs_SPLIThalf.m)
%   big_data_cond_ind: an array indicating which task condition each matrix comes from (1-5)
%   big_data_sub_ind: an array indicating which subject each matrix comes from (1-10)
%   big_data_sess_ind: an array indicating which session each matrix comes from (1-2 or 1-10)
%   atlas_params: a structure with atlas info (see atlas_parameters.m)
%   outdir_orig: where data/plots should be stored (will amake a subfolder within this)
%   title_start: the names for outputs
%
% CG


% set up folder for saving outputs
outdir = [outdir_orig 'relationship_by_node/'];
if ~exist(outdir)
    mkdir(outdir);
end

% conditions
%Ctypes = {'win_sub_btwn_cond','win_sub_and_sess','win_cond_and_sub','win_cond_btwn_sub','not_cond_or_sub'};
Ctypes = {'not_cond_or_sub','win_cond_btwn_sub','win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub'};
Ctypes_short = {'group','task','individual','indbysess','indbytask'}; %short versions for figure labels

% parcel giftis - needed for creating parcel files
watershed_L = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_L_watershedmerge_0.35_tweaked.func.gii'];
watershed_R = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'];

% loop through parcels
for n = 1:atlas_params.num_rois

    % ROI lists
    in_rois = n;
    out_rois = setdiff(1:atlas_params.num_rois,in_rois); % all other ROIs (don't take diagonal)
    
    % create a linearized version of the big_data_mat, per each ROI
    for i = 1:length(big_data_cond_ind)
        big_data_lin(i,:) = big_data_mat(in_rois,out_rois,i);
    end
    
    % take correlations for that row across different conditions
    big_data_corr = atanh(corr(big_data_lin')); 
    
    % measure similarity across conditions matched on different factors
    [simvals subvals condvals] = calc_corr_sets2(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind);
    
    % for each condition type, calculate effect (i.e., subtract out baseline)
    tot = 0;
    for ct = 1:length(Ctypes)
        allvals.(Ctypes{ct})(n) = simvals.(Ctypes{ct});
        if ct == 1 % group - no baseline
            normvals.(Ctypes{ct})(n) = simvals.(Ctypes{ct});
        elseif ct == 2 || ct == 3 %task and ind normed by group
            normvals.(Ctypes{ct})(n) = simvals.(Ctypes{ct}) - simvals.(Ctypes{1});
        else % ind&task and ind&sess normed by ind
            normvals.(Ctypes{ct})(n) = simvals.(Ctypes{ct}) - simvals.(Ctypes{3});
        end
        tot = tot + normvals.(Ctypes{ct})(n);
    end
    
    % measure normalized relative effect by finding proporition of total
    for ct = 1:length(Ctypes)
        normRelvals.(Ctypes{ct})(n) = normvals.(Ctypes{ct})(n)/tot;
    end
end

% assign data to parcel
for ct = 1:length(Ctypes)    
    outputname = [title_start '_' Ctypes_short{ct}];
    assign_data_to_parcel_cifti(allvals.(Ctypes{ct})',watershed_L,watershed_R,outdir,outputname);
end

% and a few others that I think could be interesting:
% compared with "baseline"
for ct = 1:length(Ctypes)
    assign_data_to_parcel_cifti(normvals.(Ctypes{ct})',watershed_L,watershed_R,outdir,[title_start '_' Ctypes_short{ct} '_NORM']);
end

%%% and plot as normalized relative effect sizes:
for ct = 1:length(Ctypes)
    assign_data_to_parcel_cifti(normRelvals.(Ctypes{ct})',watershed_L,watershed_R,outdir,[title_start '_' Ctypes_short{ct} '_NORMrel']);
end


end

function [simvals subvals condvals] = calc_corr_sets2(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind)

% To calc:
% avg within sub
% avg within cond
% avg within sub within cond
% avg within sub btwn cond
% avg btwn sub within cond

%if ~strcmp(distType,'euclidean')
%    big_data_corr = atanh(big_data_corr);
%end
big_data_corr(logical(eye(size(big_data_corr)))) = nan;
maskmat = ones(size(big_data_corr));
maskmat = logical(triu(maskmat,1));

conds = unique(big_data_cond_ind);
subs = unique(big_data_sub_ind);
sessions = unique(big_data_sess_ind);

cond_ind_mat = zeros(size(big_data_corr));
condspec_ind_mat = zeros(size(big_data_corr));
for c = conds
    cond_ind_mat(big_data_cond_ind == c,big_data_cond_ind==c) = 1;
    condspec_ind_mat(big_data_cond_ind == c,:) = c;
end
cond_ind_mat(logical(eye(size(big_data_corr)))) = 0;
cond_ind_mat = logical(cond_ind_mat);
condspec_ind_mat(logical(eye(size(big_data_corr)))) = 0;

% and do by sessions too... a bit more annoying
sess_ind_mat = zeros(size(big_data_corr));
for se = sessions
    sess_ind_mat(big_data_sess_ind == se,big_data_sess_ind==se) = 1;
end
sess_ind_mat(logical(eye(size(big_data_corr)))) = 0;
sess_ind_mat = logical(sess_ind_mat);

% and by subject
sub_ind_mat = zeros(size(big_data_corr));
subspec_ind_mat = zeros(size(big_data_corr));
for s = subs
    sub_ind_mat(big_data_sub_ind == s,big_data_sub_ind==s) = 1;
    subspec_ind_mat(big_data_sub_ind == s,:) = s;    
end
sub_ind_mat(logical(eye(size(big_data_corr)))) = 0;
sub_ind_mat = logical(sub_ind_mat);
subspec_ind_mat(logical(eye(size(big_data_corr)))) = 0;
%subspec_ind_mat = subspec_ind_mat.*maskmat;

% compute all of the values I want, per subject
for s=subs
    subvals.all_avg(s) = mean(nanmean(big_data_corr(big_data_sub_ind == s,:)));
    subvals.win_cond(s) =  mean(big_data_corr(logical(cond_ind_mat.*(subspec_ind_mat ==s))));
    subvals.win_sub(s) = mean(big_data_corr(logical(sub_ind_mat.*(subspec_ind_mat == s))));
    subvals.win_cond_and_sub(s) = mean(big_data_corr(logical(cond_ind_mat.*sub_ind_mat.*(subspec_ind_mat==s))));
    subvals.win_sub_and_sess(s) = mean(big_data_corr(logical(sess_ind_mat.*sub_ind_mat.*(subspec_ind_mat==s))));
    subvals.win_cond_btwn_sub(s) = mean(big_data_corr(logical(cond_ind_mat.*~sub_ind_mat.*(subspec_ind_mat==s))));
    subvals.win_sub_btwn_cond(s) = mean(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat.*(subspec_ind_mat==s))));
    subvals.win_sub_btwn_cond_and_sess(s) = mean(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat.*~sess_ind_mat.*(subspec_ind_mat==s))));
    subvals.not_cond_or_sub(s) = mean(big_data_corr(logical(~cond_ind_mat.*~sub_ind_mat.*(subspec_ind_mat==s))));
end

% compute all of the values I want, per task
for c = conds
    condvals.all_avg(c) = mean(nanmean(big_data_corr(big_data_cond_ind == c,:)));
    condvals.win_cond(c) =  mean(big_data_corr(logical(cond_ind_mat.*(condspec_ind_mat == c))));
    condvals.win_sub(c) = mean(big_data_corr(logical(sub_ind_mat.*(condspec_ind_mat == c))));
    condvals.win_cond_and_sub(c) = mean(big_data_corr(logical(cond_ind_mat.*sub_ind_mat.*(condspec_ind_mat==c))));
    condvals.win_sub_and_sess(c) = mean(big_data_corr(logical(sess_ind_mat.*sub_ind_mat.*(condspec_ind_mat==c))));
    condvals.win_cond_btwn_sub(c) = mean(big_data_corr(logical(cond_ind_mat.*~sub_ind_mat.*(condspec_ind_mat==c))));
    condvals.win_sub_btwn_cond(c) = mean(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat.*(condspec_ind_mat==c))));
    condvals.win_sub_btwn_cond_and_sess(c) = mean(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat.*~sess_ind_mat.*(condspec_ind_mat==c))));
    condvals.not_cond_or_sub(c) = mean(big_data_corr(logical(~cond_ind_mat.*~sub_ind_mat.*(condspec_ind_mat==c))));
end

% take means
types = {'all_avg','win_cond','win_sub','win_cond_and_sub','win_sub_and_sess',...
    'win_cond_btwn_sub','win_sub_btwn_cond','win_sub_btwn_cond_and_sess','not_cond_or_sub'};

for t = 1:length(types)
    simvals.(types{t}) = mean(subvals.(types{t})(subs));
    sename = [types{t} '_se'];
    simvals.(sename) = std(subvals.(types{t})(subs))/sqrt(length(subs));
end

      

end