function calc_linearmodel_byparcel(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,nameStart)
% function calc_linearmodel_byparcel(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,nameStart)
% function similar to calc_linearmodel_byedge, but does it per parcel for
% 1D maps (e.g., activation) rather than 2D matrices (e.g., FC)
%
% Inputs
%   big_data_mat: a node x type matrix, where type = a particular
%      subject-session-task grouping for a dataset (made in similarity_figs_SPLIThalf.m)
%   big_data_cond_ind: an array indicating which task condition each matrix comes from (1-5)
%   big_data_sub_ind: an array indicating which subject each matrix comes from (1-10)
%   big_data_sess_ind: an array indicating which session each matrix comes from (1-2 or 1-10)
%   atlas_params: a structure with atlas info (see atlas_parameters.m)
%   outdir: where data/plots should be stored
%   nameStart: the names for outputs
%
% CG

xVars = [big_data_sub_ind', big_data_sess_ind', big_data_cond_ind'];
mat_names = {'sub','sess','task','subXsess','subXtask','sessXtask'};

Fstat_mat = zeros(size(big_data_mat,1),6);
Pval_mat = zeros(size(big_data_mat,1),6);
VarExpBiased_mat = zeros(size(big_data_mat,1),6);
VarExp_mat = zeros(size(big_data_mat,1),6);
ModelVarExp_mat = zeros(size(big_data_mat),1);

%matlabpool open 4
for i1 = 1:size(big_data_mat,1)
    disp(['ROI: ' num2str(i1)]);
    yVars = squeeze(big_data_mat(i1,:));
    
    mdl = LinearModel.fit(xVars,yVars,'interactions','VarNames',{'sub','sess','task','FC'},'CategoricalVars',[1,2,3]); %,'RobustOpts','on');
    
    Fstat_mat(i1,:) = double(mdl.anova.F(1:end-1));
    Pval_mat(i1,:) = double(mdl.anova.pValue(1:end-1));
    
    % omega-sqared formula:
    % http://daniellakens.blogspot.com/2015/06/why-you-should-use-omega-squared.html
    %VarExp_mat(i1,i2,:) = (mdl.anova.DF(1:end-1).*(mdl.anova.MeanSq(1:end-1)-mdl.anova.MeanSq(end)))./(mdl.anova.DF(1:end-1).*mdl.anova.MeanSq(1:end-1)+((size(xVars,1)-mdl.anova.DF(1:end-1)).*mdl.anova.MeanSq(end)));
    % other formula:
    % http://www.theanalysisfactor.com/calculate-effect-size/
    %VarExp_mat(i1,i2,:) = double((mdl.anova.SumSq(1:end-1)-(mdl.anova.DF(1:end-1).*mdl.anova.MeanSq(end)))./(sum(mdl.anova.SumSq)-mdl.anova.MeanSq(end)));
    % https://egret.psychol.cam.ac.uk/statistics/local_copies_of_sources_Cardinal_and_Aitken_ANOVA/glm_effectsize.htm
    % http://www.statisticshowto.com/omega-squared/
    VarExp_mat(i1,:) = double((mdl.anova.SumSq(1:end-1)-(mdl.anova.DF(1:end-1).*mdl.anova.MeanSq(end)))./(sum(mdl.anova.SumSq)+mdl.anova.MeanSq(end)));
    
    % simpler eta-squared formula (but biased)
    VarExpBiased_mat(i1,:) = double(mdl.anova.SumSq(1:end-1))./sum(double(mdl.anova.SumSq));
    
    % variance explained for the model
    ModelVarExp_mat(i1) = 1-double(mdl.anova.SumSq(end))./sum(double(mdl.anova.SumSq));
    
end
%matlabpool close

save([outdir 'paracelwise_model_parameters_' nameStart '.mat'],'Fstat_mat','Pval_mat','VarExp_mat','ModelVarExp_mat','VarExpBiased_mat','mat_names');

watershed_L = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_L_watershedmerge_0.35_tweaked.func.gii'];
watershed_R = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'];
addpath('/data/cn/data1/scripts/CIFTI_RELATED/Resources/');

outputname = ['parcelwise_model_variance_explained_' nameStart];
assign_data_to_parcel_cifti(ModelVarExp_mat,watershed_L,watershed_R,outdir,outputname);
Batch_wb_image_capture([outdir outputname '.dtseries.nii'],[outdir outputname],'palette','videen_style','colorscaling',[0 1 0 -1]);

Fthresholds = [0 20; 0 20; 0 20; 0 20; 0 20; 0 20];
Vthresholds = [0 0.3; 0 0.3; 0 .3; 0 0.3; 0 0.3; 0 0.3];

% save out variables to cifti parcel files
for m = 1:length(mat_names)
    outputname = ['parcelwise_F_' mat_names{m} '_' nameStart];
    assign_data_to_parcel_cifti(Fstat_mat(:,m),watershed_L,watershed_R,outdir,outputname);
    Batch_wb_image_capture([outdir outputname '.dtseries.nii'],[outdir outputname],'palette','videen_style','colorscaling',[Fthresholds(m,1) Fthresholds(m,2) 0 -1]);

    outputname = ['parcelwise_Pvals_' mat_names{m} '_' nameStart];
    aa = squeeze(Pval_mat(:,m));
    [p_FDR p_FDR_masked] = FDR(aa);    
    assign_data_to_parcel_cifti(p_FDR,watershed_L,watershed_R,outdir,outputname);
    Batch_wb_image_capture([outdir outputname '.dtseries.nii'],[outdir outputname],'palette','videen_style','colorscaling',[0 0.05 0 -1]);

    outputname = ['parcelwise_Omegasquared_' mat_names{m} '_' nameStart];
    assign_data_to_parcel_cifti(VarExp_mat(:,m),watershed_L,watershed_R,outdir,outputname);
    Batch_wb_image_capture([outdir outputname '.dtseries.nii'],[outdir outputname],'palette','videen_style','colorscaling',[Vthresholds(m,1) Vthresholds(m,2) 0 -1]);

    outputname = ['parcelwise_Etasquared_' mat_names{m} '_' nameStart];
    assign_data_to_parcel_cifti(VarExpBiased_mat(:,m),watershed_L,watershed_R,outdir,outputname);
    Batch_wb_image_capture([outdir outputname '.dtseries.nii'],[outdir outputname],'palette','videen_style','colorscaling',[Vthresholds(m,1) Vthresholds(m,2) 0 -1]);

end



end