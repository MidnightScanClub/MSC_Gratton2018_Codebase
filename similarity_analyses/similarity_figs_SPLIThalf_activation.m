function similarity_figs_bysess_activation()
% function similarity_figs_bysess_activation(matchType)
% run similarity analyses to evaluate similarity across activation maps
%   this version does analyses on 5-session groups (split half)
%
% CG

matchType = 'orig'; %only this version is coded currently

% directory information
topdir = '/data/nil-bluearc/GMT/Caterina/TaskFC/FC_Parcels/';
outdir = [topdir 'similarity_analyses_' matchType '_activation/'];
if ~exist(outdir)
    mkdir(outdir);
end

% set up some variables
tasks_all = {'mixed','motor','mem'};
cond_types_orig = {{'AllGlass','AllSemantic'},{'AllCond'},{'AllCond'}};
cond_types_all = {{'AllGlass','AllSemantic'},{'AllMotor'},{'AllMem'}};
cond_types_list = {'Coherence','Semantic','Motor','Memory'};
%cond_types_all = {{'rest'},{'AllGlass','AllSemantic'},{'AllMotor'},{'Face','Scene','Word'}};
atlas_params = atlas_parameters('Parcels','/data/cn5/caterina/Atlases/Evan_parcellation/');

% colors
tmp_colors = distinguishable_colors(30);
%tmp_colors(1:3,:) = [1 0 0; 0 1 0; 0 0 1];
%color_list.conds = tmp_colors(1:length(cond_types_list),:);
color_list.conds = tmp_colors(13:length(cond_types_list)+12,:);

% set up some variables
switch matchType
    case 'orig'
    subjects = [1:7,9:10]; %1:10;
    subject_list = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
    color_list.subs = [0 0 0; 0.9 0.9 0; 0 1 0; 1 0 0; 0 0 1; 0.2 1 1; 1 0 1; 0 0.6 0.6; 1 0.5 0];
    case 'match'
        %subjects = [1:8];
end
%sessions = 1:10;
sessions = 1:2;
all_minusUnl = [2:13];
cont_net = [4:8];
proc_net = [3,9,10,11];

% load activation maps (created with create_task_activation_mats.m)
for t = 1:length(tasks_all)
    task_data.(tasks_all{t}) = load([topdir tasks_all{t} '/percent_change_allsubs.mat']);
end


count = 1; count2 = 1;
bcount = 1; bcount2 = 1;
for s = subjects
    cond_count = 1; cond_count2 = 1;
    for t = 1:length(tasks_all)
        
        %identify conds
        conds = cond_types_orig{t};
        for c = 1:length(conds)
            
            % idenfity good runs
            good_runs = [];
            for i = sessions
                badrun = sum(isnan(task_data.(tasks_all{t}).allsubs_parcel_act.(conds{c}){s}{i}));
                if badrun == 0
                    good_runs = [good_runs i];
                end
            end
            nruns = length(good_runs);
            
            % create a large matrix that is node x [type] where type
            % is an activation map from a specific subject, split-half session, and task
            big_data_mat(:,count:(count+nruns-1)) = squeeze(task_data.(tasks_all{t}).allsubs_parcel_act_concat_halves.(conds{c})(:,s,good_runs));
            big_data_sub_ind(count:(count+nruns-1)) = s;
            big_data_cond_ind(count:(count+nruns-1)) = cond_count; %t;
            big_data_sess_ind(count:(count+nruns-1)) = good_runs;
            count = count+nruns;
            cond_count = cond_count+1;
            
        end
        
    end
end

save([outdir 'data_mats.mat'],'big_data_mat','big_data_sub_ind','big_data_cond_ind');

% do linear regression model by parcel - TAKES A LITTLE WHILE
% calc_linearmodel_byparcel(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,'Raw');

% make regular similarity figures, organized by subject
[h1 h2 simvals.raw subvals.raw condvals.raw] = make_corr_relationship_SUBSET_figs(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,all_minusUnl,atlas_params,'AllminusUnl',find(diff(big_data_sub_ind)));
save_fig(h1,[outdir 'AllminusUnl_raw_raster.pdf']);
save_fig(h2,[outdir 'AllminusUnl_raw.pdf']);

% do the same thing for control networks only
[h1 h2 simvals.raw_ctrl subvals.raw_ctrl condvals.raw_ctrl] = make_corr_relationship_SUBSET_figs(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,cont_net,atlas_params,'ControlOnly',find(diff(big_data_sub_ind)));
save_fig(h1,[outdir 'ControlOnly_raw_raster.pdf']);
save_fig(h2,[outdir 'ControlOnly_raw.pdf']);

% and for processing networks only
[h1 h2 simvals.raw_proc subvals.raw_proc condvals.raw_proc] = make_corr_relationship_SUBSET_figs(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,proc_net,atlas_params,'ProcOnly',find(diff(big_data_sub_ind)));
save_fig(h1,[outdir 'ProcOnly_raw_raster.pdf']);
save_fig(h2,[outdir 'ProcOnly3_raw.pdf']);

% save out results
save([outdir 'similarity_vals.mat'],'subvals','subjects','condvals','simvals');

% make summary figurse
similarity_quant_fig3(subvals,matchType,'corr',subjects,0);
save_fig(gcf,[outdir 'Quant_raw_diff_FULLonly.pdf']);
similarity_quant_fig2(subvals,matchType,'corr',subjects,0);
save_fig(gcf,[outdir 'Quant_raw_diff_FULLonly_OLDver.pdf']);
similarity_quant_fig2(subvals,matchType,'corr',subjects,1);
save_fig(gcf,[outdir 'Quant_raw_diff.pdf']);
%similarity_quant_fig(simvals,matchType)
%save_fig(gcf,[outdir 'Quant_raw_diff.pdf']);
similarity_quant_fig_bysub(subvals,matchType,subjects,subject_list,'corr',color_list);
save_fig(gcf,[outdir 'Quant_raw_diff_bysub.pdf']);
similarity_quant_fig_bycond(condvals,matchType,cond_types_all,'corr',color_list);
save_fig(gcf,[outdir 'Quant_raw_diff_bycond.pdf']);
close('all');

end

function make_corr_relationship_figs(big_data_mat,big_data_cond_ind,big_data_sub_ind)

clear big_data_lin;
maskmat = ones(size(big_data_mat,1));
maskmat = triu(maskmat,1);
for i = 1:length(big_data_cond_ind)
    temp = big_data_mat(:,:,i);
    big_data_lin(i,:) = temp(logical(maskmat));
end
big_data_corr = corr(big_data_lin');

figure('Position',[1 1 1600 400]);
subplot(1,5,1:3);
imagesc(big_data_lin,[-1 1]);
xlabel('edges'); ylabel('subs,conds'); title('linearized corrmats');
subplot(1,5,4);
imagesc(big_data_cond_ind');
title('task');
subplot(1,5,5);
imagesc(big_data_sub_ind');
title('subject');


figure('Position',[1 1 1400 500]);
subplot(1,4,1:2);
imagesc(big_data_corr,[-0.4 1]);
xlabel('all subs & conds'); ylabel('all subs & conds'); title('correlation between matrices');
colorbar();
subplot(1,4,3);
imagesc(big_data_cond_ind');
title('task');
subplot(1,4,4);
imagesc(big_data_sub_ind');
title('subject');

end

function [h1 h2 simvals subvals condvals] = make_corr_relationship_SUBSET_figs(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,net_subset,atlas_params,fulltitle,addlines,varargin)

net_inds = [];
for n = 1:length(net_subset)
    net_inds = [net_inds atlas_params.mods{net_subset(n)}];
end

%maskmat = ones(length(net_inds));
%maskmat = triu(maskmat,1);
for i = 1:length(big_data_cond_ind)
    %temp = big_data_mat(net_inds,net_inds,i);
    big_data_lin(i,:) = big_data_mat(net_inds,i); %temp(logical(maskmat));
end
big_data_corr = corr(big_data_lin');
%simvals = calc_corr_sets(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind);
[simvals subvals condvals] = calc_corr_sets2(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind);

h1 = figure('Position',[1 1 1600 400]);
subplot(1,5,1:3);
imagesc(big_data_lin,[-1 1]);
xlabel('parcels'); ylabel('subs,conds'); title('linearized corrmats');
hline_new(addlines+0.5,'k',1);
axis square;
subplot(1,5,4);
imagesc(big_data_cond_ind');
title('task');
subplot(1,5,5);
imagesc(big_data_sub_ind');
title('subject');
figtitle(fulltitle);

h2 = figure('Position',[1 1 1400 500]);
subplot(1,4,1:2);
if isempty(varargin)
    imagesc(big_data_corr,[0 1.5]);
else
    imagesc(big_data_corr,[0 0.75]);
end
xlabel('all subs & conds'); ylabel('all subs & conds'); title('correlation between matrices');
hline_new(addlines+0.5,'k',1);
vline_new(addlines+0.5,'k',1);
colorbar();
axis square;
subplot(1,4,3);
imagesc(big_data_cond_ind');
title('task');
subplot(1,4,4);
imagesc(big_data_sub_ind');
title('subject');
figtitle(fulltitle);
end

function [simvals] = corr_relationship_bynod(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,title_start)

Ctypes = {'win_sub_btwn_cond','win_sub_and_sess','win_cond_and_sub','win_cond_btwn_sub','not_cond_or_sub'};
watershed_L = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_L_watershedmerge_0.35_tweaked.func.gii'];
watershed_R = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'];

for n = 1:atlas_params.num_rois
%    in_rois = n;
%    out_rois = setdiff(1:atlas_params.num_rois,in_rois);
%    
%for i = 1:length(big_data_cond_ind)
%    big_data_lin(i,:) = big_data_mat(in_rois,out_rois,i);
%end
big_data_lin = big_data_mat';

big_data_corr = corr(big_data_lin');
simvals = calc_corr_sets(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind);
for ct = 1:length(Ctypes)
    allvals.(Ctypes{ct})(n) = simvals.(Ctypes{ct});
end
end

for ct = 1:length(Ctypes)    
    outputname = [title_start '_' Ctypes{ct}];
    assign_data_to_parcel_cifti(allvals.(Ctypes{ct})',watershed_L,watershed_R,outdir,outputname);
end

% and a few others that I think could be interesting:
% compared with "baseline"
assign_data_to_parcel_cifti(allvals.win_sub_btwn_cond' - allvals.not_cond_or_sub',watershed_L,watershed_R,outdir,[title_start '_subOnly_vs_none']);
assign_data_to_parcel_cifti(allvals.win_sub_and_sess' - allvals.not_cond_or_sub',watershed_L,watershed_R,outdir,[title_start '_subANDsess_vs_none']);
assign_data_to_parcel_cifti(allvals.win_cond_and_sub' - allvals.not_cond_or_sub',watershed_L,watershed_R,outdir,[title_start '_subANDcond_vs_none']);
assign_data_to_parcel_cifti(allvals.win_cond_btwn_sub' - allvals.not_cond_or_sub',watershed_L,watershed_R,outdir,[title_start '_condOnly_vs_none']);

% sub&cond - subOnly; sub&sess-subOnly
assign_data_to_parcel_cifti(allvals.win_cond_and_sub' - allvals.win_sub_btwn_cond',watershed_L,watershed_R,outdir,[title_start '_subANDcond_vs_subOnly']);
assign_data_to_parcel_cifti(allvals.win_sub_and_sess' - allvals.win_sub_btwn_cond',watershed_L,watershed_R,outdir,[title_start '_subANDsess_vs_subOnly']);

% sub vs. cond
assign_data_to_parcel_cifti(allvals.win_sub_btwn_cond' - allvals.win_cond_btwn_sub',watershed_L,watershed_R,outdir,[title_start '_subOnly_vs_condOnly']);


end

function similarity_quant_fig(simvals,matchType)

types = {'raw','raw_ctrl','raw_proc'};
type_names = {'full','ctrl','proc'};
comps = {'win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub','win_cond_btwn_sub','not_cond_or_sub'};
comps_se = {'win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se','win_cond_btwn_sub_se','not_cond_or_sub_se'};
comps_short = {'subOnly','sub&sess','sub&cond','condOnly','~condOrSub'};

colors = distinguishable_colors(length(types));

figure('Position',[1 1 800 800]);
subplot(1,2,1)
for t = 1:length(types)
    plotvals = []; plotvals_se = [];
    for c = 1:length(comps)
        plotvals = [plotvals simvals.(types{t}).(comps{c})];
        plotvals_se = [plotvals_se simvals.(types{t}).(comps_se{c})];
    end
    %plot(1:length(plotvals),plotvals,'.','Color',colors(t,:),'MarkerSize',15); hold on;
    errorbar(1:length(plotvals),plotvals,plotvals_se,'.','Color',colors(t,:),'MarkerSize',5); hold on;
end
ylabel('r (z)');
xlim([0.5 length(plotvals)+0.5]);
switch matchType
    case 'orig'
        ylim([0 2]); %ylim([0.5 1.6]);
    case 'match'
        ylim([0.4 1.2]);
end
legend(type_names);
axis square;
set(gca,'XTick',1:length(plotvals),'Xticklabels',comps_short);
xticklabel_rotate;
title('Similarity: Task conds vs. Baseline');



end

function simvals = calc_corr_sets(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind)

% To calc:
% avg within sub
% avg within cond
% avg within sub within cond
% avg within sub btwn cond
% avg btwn sub within cond

big_data_corr = atanh(big_data_corr);
%big_data_corr(logical(eye(size(big_data_corr)))) = nan;
maskmat = ones(size(big_data_corr));
maskmat = logical(triu(maskmat,1));

conds = unique(big_data_cond_ind);
subs = unique(big_data_sub_ind);
sessions = unique(big_data_sess_ind);

cond_ind_mat = zeros(size(big_data_corr));
for c = conds
    cond_ind_mat(big_data_cond_ind == c,big_data_cond_ind==c) = 1;
end
cond_ind_mat(logical(eye(size(big_data_corr)))) = 0;
cond_ind_mat = logical(cond_ind_mat.*maskmat);
    
sub_ind_mat = zeros(size(big_data_corr));
for s = subs
    sub_ind_mat(big_data_sub_ind == s,big_data_sub_ind==s) = 1;
end
sub_ind_mat(logical(eye(size(big_data_corr)))) = 0;
sub_ind_mat = logical(sub_ind_mat.*maskmat);

% and do by sessions too... a bit more annoying
sess_ind_mat = zeros(size(big_data_corr));
for se = sessions
    sess_ind_mat(big_data_sess_ind == se,big_data_sess_ind==se) = 1;
end
sess_ind_mat(logical(eye(size(big_data_corr)))) = 0;
sess_ind_mat = logical(sess_ind_mat.*maskmat);


simvals.all_avg = mean(big_data_corr(maskmat));
simvals.win_cond = mean(big_data_corr(cond_ind_mat));
simvals.win_sub = mean(big_data_corr(sub_ind_mat));
simvals.win_cond_and_sub = mean(big_data_corr(logical(cond_ind_mat.*sub_ind_mat)));
simvals.win_sub_and_sess = mean(big_data_corr(logical(sess_ind_mat.*sub_ind_mat)));
simvals.win_cond_btwn_sub = mean(big_data_corr(logical(cond_ind_mat.*~sub_ind_mat)));
simvals.win_sub_btwn_cond = mean(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat)));
simvals.not_cond_or_sub = mean(big_data_corr(logical(~cond_ind_mat.*~sub_ind_mat.*maskmat)));

% and compute error bars (?? these are a bit weird, since not really
% independent
simvals.all_avg_se = std(big_data_corr(maskmat))/sqrt(sum(sum(maskmat)));
simvals.win_cond_se = std(big_data_corr(cond_ind_mat))/sqrt(sum(sum(cond_ind_mat)));
simvals.win_sub_se = std(big_data_corr(sub_ind_mat))/sqrt(sum(sum(sub_ind_mat)));
simvals.win_cond_and_sub_se = std(big_data_corr(logical(cond_ind_mat.*sub_ind_mat)))/sqrt(sum(sum(cond_ind_mat.*sub_ind_mat)));
simvals.win_sub_and_sess_se = std(big_data_corr(logical(sess_ind_mat.*sub_ind_mat)))/sqrt(sum(sum(sess_ind_mat.*sub_ind_mat)));
simvals.win_cond_btwn_sub_se = std(big_data_corr(logical(cond_ind_mat.*~sub_ind_mat)))/sqrt(sum(sum(cond_ind_mat.*~sub_ind_mat)));
simvals.win_sub_btwn_cond_se = std(big_data_corr(logical(sub_ind_mat.*~cond_ind_mat)))/sqrt(sum(sum(~cond_ind_mat.*sub_ind_mat)));
simvals.not_cond_or_sub_se = std(big_data_corr(logical(~cond_ind_mat.*~sub_ind_mat.*maskmat)))/sqrt(sum(sum(~cond_ind_mat.*~sub_ind_mat.*maskmat)));


end

function [simvals subvals condvals] = calc_corr_sets2(big_data_corr,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind)

% To calc:
% avg within sub
% avg within cond
% avg within sub within cond
% avg within sub btwn cond
% avg btwn sub within cond

big_data_corr = atanh(big_data_corr);
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

function similarity_quant_fig_bycond(condvals,matchType,condlist,distType,color_list)

% convert subjects to cell array
count = 0;
for c1 = 1:length(condlist)
    condsets = condlist{c1};
    for c2 = 1:length(condsets)
        count = count + 1;
        conditions{count} = condsets{c2};
    end
end

types = {'raw'};
type_names = {'full'};
comps = {'not_cond_or_sub','win_cond_btwn_sub','win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub'};
comps_se = {'not_cond_or_sub_se','win_cond_btwn_sub_se','win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se'};
comps_short = {'group','task','indiv','indiv*sess','indiv*task'};

%colors = distinguishable_colors(length(conditions));
colors = color_list.conds;

figure('Position',[1 1 800 800]);

subplot(1,2,1)
for s = 1:count
for t = 1:length(types)
    plotvals = [];
    for c = 1:length(comps)
        plotvals = [plotvals condvals.(types{t}).(comps{c})(s)];
    end
    plot(1:length(plotvals),plotvals,'.','Color',colors(s,:),'MarkerSize',15); hold on;
end
end
ylabel(distType);
xlim([0.5 length(plotvals)+0.5]);
%switch matchType
%    case 'orig'
%        ylim([0 2]); %ylim([0.5 1.6]);
%    case 'match'
%        ylim([0.4 1.2]);
%end
legend(conditions);
axis square;
set(gca,'XTick',1:length(plotvals),'Xticklabels',comps_short);
xticklabel_rotate;
title('Similarity: Raw Matrix');


end


function similarity_quant_fig_bysub(subvals,matchType,sublist,subject_list,distType,color_list)

% convert subjects to cell array
%for s = 1:length(sublist)
%    subjects{s} = sprintf('MSC%02d',sublist(s));
%end

types = {'raw'};
type_names = {'full'};
comps = {'not_cond_or_sub','win_cond_btwn_sub','win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub'};
comps_se = {'not_cond_or_sub_se','win_cond_btwn_sub_se','win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se'};
comps_short = {'group','task','indiv','indiv*sess','indiv*task'};

%colors = distinguishable_colors(max(sublist));
colors = color_list.subs;

figure('Position',[1 1 800 800]);

subplot(1,2,1)
for s = 1:length(sublist)
for t = 1:length(types)
    plotvals = [];
    for c = 1:length(comps)
        plotvals = [plotvals subvals.(types{t}).(comps{c})(sublist(s))];
    end
    plot(1:length(plotvals),plotvals,'.','Color',colors(s,:),'MarkerSize',15); hold on;
end
end
ylabel(distType);
xlim([0.5 length(plotvals)+0.5]);
%switch matchType
%    case 'orig'
%        ylim([0 2]); %ylim([0.5 1.6]);
%    case 'match'
%        ylim([0.4 1.2]);
%end
legend(subject_list); %legend(subjects);
set(gca,'XTick',1:length(plotvals),'Xticklabels',comps_short);
xticklabel_rotate;
axis square;
title('Similarity: Raw Matrix');


end

function make_boxplot(subvals,comps,types,subjects,distType,matchType,comps_short,type_names)

%colors = distinguishable_colors(length(types));
colors = [0 0 0; 1 0 0; 0 0 1; 0 1 0];

plotvals = []; groups = []; positions = []; colors_all = [];
count = 1; start_pos = 1;
for c = 1:length(comps)
    for t = 1:length(types)
        plotvals = [plotvals; subvals.(types{t}).(comps{c})(subjects)];
        groups = [groups count];
        positions = [positions start_pos+0.2*(t-1)];
        count = count +1;
        colors_all = [colors_all; colors(t,:)];
    end
    start_pos = start_pos + 1;
end
%bh = boxplot(plotvals','colors',colors(t,:),'symbol','k+','positions',[1:size(plotvals,1)]+0.25*(t-1)); hold on; %'boxstyle','filled'); hold on;
bh = boxplot(plotvals',groups,'positions',positions,'symbol','k+','colors',colors_all); hold on; %'boxstyle','filled'); hold on;
%set(bh,'linewidth',1.5/length(types));
%h = findobj(gca,'Tag','Box');
%for j = 1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),colors_all(end-j+1,:),'FaceAlpha',0.25*length(types));
%end
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-');
set(findobj(gca,'type','line'),'linew',3);

switch distType
    case 'corr'
        ylabel(distType);
        lim1 = [0 2.5]; %[0 2];
        lim2 = [0.4 1.2];
    case '1minuscorr'
        ylabel(distType);
        lim1 = [0 1];
        lim2 = [0 1];
    case 'euclidean'
        ylabel(distType);
        lim1=[0 50];
        lim2=[0 50];
    case 'ratioComp'
        ylabel('Proportion Removed: (raw - diff)/raw')
        lim1 = [0 1];
        lim2 = [0 1];
end

switch matchType
    case 'orig'
        ylim(lim1); %ylim([0.5 1.6]);
    case 'match'
        ylim(lim2);
end
legend(type_names);
%axis square;
ax = axis;
%axis(axis);
Yl = ax(3:4); 
xlim([positions(1)-0.5 positions(end)+0.5]);
if length(types) > 1
    Xt = positions([1:length(comps)]*length(types)-1);
else
    Xt = positions([1:length(comps)]*length(types));
end    
set(gca,'XTick',Xt,'XTickLabel',{' '});
t = text(Xt,Yl(1)*ones(1,length(Xt)),comps_short);
set(t,'HorizontalAlignment','right','VerticalAlignment','middle','Rotation',90);
%xticklabel_rotate;
end

function similarity_quant_fig2(subvals,matchType,distType,subjects,doAllTypes)

if doAllTypes
    types = {'raw','raw_ctrl','raw_proc'};
    type_names = {'full','ctrl','proc'};
else
    types = {'raw'};
    type_names = {'full'};
end
%comps = {'win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub','win_cond_btwn_sub','not_cond_or_sub'};
%comps_se = {'win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se','win_cond_btwn_sub_se','not_cond_or_sub_se'};
%comps_short = {'subOnly','sub&sess','sub&cond','condOnly','btwn'};
comps = {'not_cond_or_sub','win_cond_btwn_sub','win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub'};
comps_se = {'not_cond_or_sub_se','win_cond_btwn_sub_se','win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se'};
comps_short = {'group','task','indiv','indiv*sess','indiv*task'};


%figure('Position',[1 1 800 400]);
figure('Position',[1 1 400+100*length(types) 900]);
subplot(2,3,1:2)
make_boxplot(subvals,comps,types,subjects,distType,matchType,comps_short,type_names);
title('Similarity: Raw Matrix');

% stacked bar plot for raw, relative and normalized
subplot(2,3,3); 
plotvals = [];
for t = 1:length(types)
    for c = 1:length(comps)
        if c == 1 % group needs no normalization
            plotvals(c,t) = mean(subvals.(types{t}).(comps{c})(subjects));
        elseif c == 2 || c == 3 % task and ind normalized relative to group
            plotvals(c,t) = mean(subvals.(types{t}).(comps{c})(subjects)) - mean(subvals.(types{t}).(comps{1})(subjects));
        else % ind*ses and ind*task normalized relative to ind
            plotvals(c,t) = mean(subvals.(types{t}).(comps{c})(subjects)) - mean(subvals.(types{t}).(comps{3})(subjects)); 
        end
    end   
    sumtot = sum(plotvals(:,t));
    plotvals_rel(:,t) = plotvals(:,t)./sumtot;
end
plotvals_rel2 = cumsum(plotvals_rel);

%plotvals_stack = [plotvals(1); diff(plotvals)];
cmap(1,:,:) = gray(size(plotvals,1));
cmap(2,:,:) = autumn(size(plotvals,1)); %cmap(2,:,2:3) = 0;
cmap(3,:,:) = winter(size(plotvals,1)); %cmap(3,:,1:2) = 0;
for c = 1:length(comps)
    for t = 1:length(types)
        bar(t,plotvals_rel2(end-c+1,t),'FaceColor',cmap(t,c,:),'EdgeColor',[0 0 0],'LineWidth',2); hold on;
    end
end
xlim([0 1+length(types)]);
set(gca,'XTick',1:length(types),'XTickLabels',types);
xticklabel_rotate();
ylim([0 1]);
%ylabel('% maximal similarity');

% do some stats - just on full for now btwn groups, and within groups
% between types
if doAllTypes
    disp('Raw, Btwn Nets')
    for t1 = 1:length(types)
        for t2 = 1:length(types)
            if t1 ~= t2
                for c1 = 1:length(comps)
                    [h p ci stats]=ttest(subvals.(types{t1}).(comps{c1})(subjects),subvals.(types{t2}).(comps{c1})(subjects));
                    disp(sprintf('%s, %s vs. %s: t(%d)=%.03f, p=%.04f',comps{c1},types{t1},types{t2},stats.df,stats.tstat,p));
                end
            end
        end
    end  
else
    disp('Raw, Full Matrix')
    for c1 = 1:length(comps)
        for c2 = 1:length(comps)
            if c1 ~= c2
                [h p ci stats]=ttest(subvals.raw.(comps{c1})(subjects),subvals.raw.(comps{c2})(subjects));
                disp(sprintf('%s vs. %s: t(%d)=%.03f, p=%.04f',comps{c1},comps{c2},stats.df,stats.tstat,p));
            end
        end
    end
end




subplot(2,3,4:5)

% stacked bar plot for raw
subplot(2,3,6); 





end

function similarity_quant_fig3(subvals,matchType,distType,subjects,doAllTypes)

if doAllTypes
    types = {'raw','raw_ctrl','raw_proc'};
    type_names = {'full','ctrl','proc'};
else
    types = {'raw'};
    type_names = {'full'};
end
%comps = {'win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub','win_cond_btwn_sub','not_cond_or_sub'};
%comps_se = {'win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se','win_cond_btwn_sub_se','not_cond_or_sub_se'};
%comps_short = {'subOnly','sub&sess','sub&cond','condOnly','btwn'};
comps = {'not_cond_or_sub','win_cond_btwn_sub','win_sub_btwn_cond_and_sess','win_sub_and_sess','win_cond_and_sub'};
comps_se = {'not_cond_or_sub_se','win_cond_btwn_sub_se','win_sub_btwn_cond_and_sess_se','win_sub_and_sess_se','win_cond_and_sub_se'};
comps_short = {'group','task','indiv','indiv*sess','indiv*task'};


%figure('Position',[1 1 800 400]);
figure('Position',[1 1 400+100*length(types) 900]);
subplot(2,3,1:2)
make_barplot(subvals,comps,types,subjects,distType,matchType,comps_short,type_names);
title('Similarity: Raw Matrix');

% do some stats - just on full for now btwn groups, and within groups
% between types
if doAllTypes
    disp('Raw, Btwn Nets')
    for t1 = 1:length(types)
        for t2 = 1:length(types)
            if t1 ~= t2
                for c1 = 1:length(comps)
                    [h p ci stats]=ttest(subvals.(types{t1}).(comps{c1})(subjects),subvals.(types{t2}).(comps{c1})(subjects));
                    disp(sprintf('%s, %s vs. %s: t(%d)=%.03f, p=%.04f',comps{c1},types{t1},types{t2},stats.df,stats.tstat,p));
                end
            end
        end
    end  
else
    disp('Raw, Full Matrix')
    for c1 = 1:length(comps)
        for c2 = 1:length(comps)
            if c1 ~= c2
                [h p ci stats]=ttest(subvals.raw.(comps{c1})(subjects),subvals.raw.(comps{c2})(subjects));
                disp(sprintf('%s vs. %s: t(%d)=%.03f, p=%.04f',comps{c1},comps{c2},stats.df,stats.tstat,p));
            end
        end
    end
end


if doAllTypes
    types = {'diff','diff_ctrl','diff_proc'};
    type_names = {'full','ctrl','proc'};
    newTypes = {'diffRatio_full','diffRatio_ctrl','diffRatio_proc'};
else
    types = {'diff'};
    type_names = {'full'};
    newTypes = {'diffRatio_full'};
end

subplot(2,3,4:5)
%make_barplot(subvals,comps,types,subjects,distType,matchType,comps_short,type_names);
%title('Similarity: Minus Rest');






end

function make_barplot(subvals,comps,types,subjects,distType,matchType,comps_short,type_names)

%colors = distinguishable_colors(length(types));
colors = [0 0 0; 1 0 0; 0 0 1; 0 1 0];

hold on;
for t = 1:length(types)
    for c = 1:length(comps)
        plotvals(c,t) = mean(subvals.(types{t}).(comps{c})(subjects));
        plotvals_se(c,t) = std(subvals.(types{t}).(comps{c})(subjects))/sqrt(length(subjects));
    end
    positions = [1:length(comps)] + (t * 0.5) - 0.5;
    h = errorbar(positions,plotvals(:,t),plotvals_se(:,t),'.','Color',colors(t,:));
    bar(positions,plotvals(:,t),'FaceColor',colors(t,:));    
    errorbar_tick(h,20);
end

switch distType
    case 'corr'
        ylabel(distType);
        lim1 = [0 1.8]; %[0 2];
        lim2 = [0 1.4];
    case '1minuscorr'
        ylabel(distType);
        lim1 = [0 1];
        lim2 = [0 1];
    case 'euclidean'
        ylabel(distType);
        lim1=[0 50];
        lim2=[0 50];
    case 'ratioComp'
        ylabel('Proportion Removed: (raw - diff)/raw')
        lim1 = [0 1];
        lim2 = [0 1];
end

switch matchType
    case 'orig'
        ylim(lim1); %ylim([0.5 1.6]);
    case 'match'
        ylim(lim2);
end
if length(types) > 1
    legend(type_names);
end
ax = axis;
Yl = ax(3:4); 
xlim([positions(1)-0.5 positions(end)+0.5]);
if length(types) > 1
    Xt = positions([1:length(comps)]*length(types)-1);
else
    Xt = positions([1:length(comps)]*length(types));
end    
set(gca,'XTick',Xt,'XTickLabel',{' '});
t = text(Xt,Yl(1)*ones(1,length(Xt)),comps_short);
set(t,'HorizontalAlignment','right','VerticalAlignment','middle','Rotation',90);
%xticklabel_rotate;

end