function calc_linearmodel_byedge(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,nameStart,varargin)
% function calc_linearmodel_byedge(big_data_mat,big_data_cond_ind,big_data_sub_ind,big_data_sess_ind,atlas_params,outdir,nameStart,varargin)
% this function does an ANOVA/linear model for each edge in a matrix
% Inputs
%   big_data_mat: a node x node x type matrix, where type = a particular
%      subject-session-task grouping for a dataset (made in similarity_figs_SPLIThalf.m)
%   big_data_cond_ind: an array indicating which task condition each matrix comes from (1-5)
%   big_data_sub_ind: an array indicating which subject each matrix comes from (1-10)
%   big_data_sess_ind: an array indicating which session each matrix comes from (1-2 or 1-10)
%   atlas_params: a structure with atlas info (see atlas_parameters.m)
%   outdir: where data/plots should be stored
%   nameStart: the names for outputs
%   varargin: mean FD array associated with each matrix
%
% CG


if nargin > 7 % include mean FD as a variable
    big_data_FD = varargin{1};
    
    big_data_sess_ind_CORR = big_data_sess_ind; % remove "extra" make up sessions
    big_data_sess_ind_CORR(big_data_sess_ind > 10) = nan;

    % a list of all the variables
    %xVars = [big_data_sub_ind', big_data_sess_ind_CORR', big_data_cond_ind', big_data_FD'];
    xVars2 = {big_data_sub_ind', big_data_sess_ind_CORR', big_data_cond_ind', big_data_FD'};

    % names for variables and their interactions
    %mat_names = {'sub','sess','task','FD','subXsess','subXtask','subXFD','sessXtask','sessXFD','taskXFD'};
    mat_names = {'sub','sess','task','FD','subXsess','subXtask','sessXtask'};
    nameStart = [nameStart '_FDincl'];
    varNames = {'sub','sess','task','FD'};
    varNames_orig = {'sub','sess','task','FD','FC'}; % for testing
    
    % terms to run
    termsNew = [1     0     0     0;...
        0     1     0     0;...
        0     0     1     0;...
        0     0     0     1;...
        1     1     0     0;...
        1     0     1     0;...
        0     1     1     0];
else
    % run this if FD is not included
    
    % fix session index to remove "extra" sessions from make up sessions
    big_data_sess_ind_CORR = big_data_sess_ind; % remove "extra" make up sessions
    big_data_sess_ind_CORR(big_data_sess_ind > 10) = nan;

    % names for variables and their interactions
    %xVars = [big_data_sub_ind', big_data_sess_ind_CORR', big_data_cond_ind'];
    xVars2 = {big_data_sub_ind', big_data_sess_ind_CORR', big_data_cond_ind'};
    mat_names = {'sub','sess','task','subXsess','subXtask','sessXtask'};
    varNames = {'sub','sess','task'};
end

% initialize some variables
Fstat_mat = zeros(size(big_data_mat,1),size(big_data_mat,2),length(mat_names));
Pval_mat = zeros(size(big_data_mat,1),size(big_data_mat,2),length(mat_names));
VarExpBiased_mat = zeros(size(big_data_mat,1),size(big_data_mat,2),length(mat_names));
VarExp_mat = zeros(size(big_data_mat,1),size(big_data_mat,2),length(mat_names));
ModelVarExp_mat = zeros(size(big_data_mat,1),size(big_data_mat,2),length(mat_names));

% loop through each edge of the matrix
for i1 = 1:size(big_data_mat,1)
    disp(['ROI: ' num2str(i1)]);
    for i2 = 1:size(big_data_mat,2)
        if i1 ~= i2
            
            % this is the edge, across all matrices
            yVars = squeeze(big_data_mat(i1,i2,:));
            
            %% Linear Model Approach%%
            %mdl = LinearModel.fit(xVars,yVars,'interactions','VarNames',{'sub','sess','task','FC'},'CategoricalVars',[1,2,3]); %,'RobustOpts','on');            
            %Fstat_mat(i1,i2,:) = double(mdl.anova.F(1:end-1));
            %Pval_mat(i1,i2,:) = double(mdl.anova.pValue(1:end-1));

            % omega-sqared formula:
            % http://www.theanalysisfactor.com/calculate-effect-size/
            % https://egret.psychol.cam.ac.uk/statistics/local_copies_of_sources_Cardinal_and_Aitken_ANOVA/glm_effectsize.htm
            % http://www.statisticshowto.com/omega-squared/
            % VarExp_mat(i1,i2,:) = double((mdl.anova.SumSq(1:end-1)-(mdl.anova.DF(1:end-1).*mdl.anova.MeanSq(end)))./(sum(mdl.anova.SumSq)+mdl.anova.MeanSq(end)));
            
            % simpler eta-squared formula (but biased)
            %VarExpBiased_mat(i1,i2,:) = double(mdl.anova.SumSq(1:end-1))./sum(double(mdl.anova.SumSq));

            % variance explained for the model
            %ModelVarExp_mat(i1,i2) = 1-double(mdl.anova.SumSq(end))./sum(double(mdl.anova.SumSq));
            
            %% ANOVAN approach: %%
            if nargin > 7 % FD included as a term
                %[p tbl stats terms] = anovan(yVars,xVars2,'varNames',varNames,'model','interaction','random',[1],'continuous',4,'display','off');
                %mdl = LinearModel.fit(xVars,yVars,'interactions','VarNames',varNames_orig,'CategoricalVars',[1,2,3]);

                % FOR TESTING
                %[p2 tbl2 stats2 terms2] = anovan(yVars,xVars2(1:3),'varNames',varNames(1:3),'model','interaction','random',[1],'display','off');
                %mdl2 = LinearModel.fit(xVars(:,1:3),yVars,'interactions','VarNames',varNames_orig([1,2,3,5]),'CategoricalVars',[1,2,3]);

                [p tbl stats terms] = anovan(yVars,xVars2,'varNames',varNames,'model', termsNew,'continuous',4,'display','off');
            else
                [p tbl stats terms] = anovan(yVars,xVars2,'varNames',varNames,'model','interaction','random',[1],'display','off');
            end
            Fstat_mat(i1,i2,:) = cell2mat(tbl(2:(length(mat_names)+1),6)); %corresponds to F column, 3 main effects and 3 interactions
            Pval_mat(i1,i2,:) = p;
            
            % omega-squared (see above for formula)
            %VarExp_mat(i1,i2,:) = double((mdl.anova.SumSq(1:end-1)-(mdl.anova.DF(1:end-1).*mdl.anova.MeanSq(end)))./(sum(mdl.anova.SumSq)+mdl.anova.MeanSq(end)));
            SumSq = cell2mat(tbl(2:end,2)); %1-6: effects, 7=error, 8=total
            DF = cell2mat(tbl(2:end,3));
            MeanSq = cell2mat(tbl(2:end,5));
            VarExp_mat(i1,i2,:) = (SumSq(1:length(mat_names))-(DF(1:length(mat_names)).*MeanSq(end-1)))./(SumSq(end)+MeanSq(end-1));
            
            % simpler eta-squared formula (but biased)
            %VarExpBiased_mat(i1,i2,:) = double(mdl.anova.SumSq(1:end-1))./sum(double(mdl.anova.SumSq));
            VarExpBiased_mat(i1,i2,:) = SumSq(1:length(mat_names))./SumSq(end);

            % variance explained for the model
            %ModelVarExp_mat(i1,i2) = 1-double(mdl.anova.SumSq(end))./sum(double(mdl.anova.SumSq));
            ModelVarExp_mat(i1,i2) = 1-SumSq(end-1)./SumSq(end);             

        end
    end
end

% save out variables
save([outdir 'edgewise_model_parameters_' nameStart '.mat'],'Fstat_mat','Pval_mat','VarExp_mat','ModelVarExp_mat','mat_names');
        
% make plots of the variance explained for the matirx
figure_corrmat_network_generic(ModelVarExp_mat,atlas_params,0,1);
colormap('jet');
title('Variance Explained, Full model');
save_fig(gcf,[outdir 'Edgewise_variance_explained_' nameStart '.pdf']);

% upper tri mask
maskmat = ones(333);
maskmat = logical(triu(maskmat,1));

% plots of the F-statistics
thresholds = [0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10; 0 10];
for m = 1:length(mat_names)
    aa = squeeze(Fstat_mat(:,:,m));
    figure_corrmat_network_generic(aa,atlas_params,thresholds(m,1),thresholds(m,2));
    colormap('jet');
    title(mat_names{m});
    save_fig(gcf,[outdir 'Edgewise_F_' mat_names{m} '_' nameStart '.pdf']);
end
close('all');

% plots of the FDR-correct p-values
for m = 1:length(mat_names)
    aa = squeeze(Pval_mat(:,:,m));
    bb = aa(maskmat);
    [p_FDR p_FDR_masked] = FDR(bb); %FDR function can be found on mathworks
    p_FDR_mat = zeros(size(aa));
    p_FDR_mat(maskmat) = p_FDR;
    p_FDR_mat = p_FDR_mat + p_FDR_mat';
    figure_corrmat_network_generic(p_FDR_mat,atlas_params,0,0.05);
    colormap(flipud(colormap('jet')));
    title(mat_names{m});
    save_fig(gcf,[outdir 'Edgewise_Pvals_' mat_names{m} '_' nameStart '.pdf']);
end
close('all');

% plots of the variance explained per variable (omega squared)
thresholds = [0 0.6; 0 0.2; 0 .2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2];
for m = 1:length(mat_names)
    aa = squeeze(VarExp_mat(:,:,m));
    figure_corrmat_network_generic(aa,atlas_params,thresholds(m,1),thresholds(m,2));
    colormap('jet');
    title(mat_names{m});
    save_fig(gcf,[outdir 'Edgewise_Omegasquared_' mat_names{m} '_' nameStart '.pdf']);
end
close('all');

% plots of the variance explained per variable (eta squared)
thresholds = [0 0.6; 0 0.2; 0 .2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2; 0 0.2];
for m = 1:length(mat_names)
    aa = squeeze(VarExpBiased_mat(:,:,m));
    figure_corrmat_network_generic(aa,atlas_params,thresholds(m,1),thresholds(m,2));
    colormap('jet');
    title(mat_names{m});
    save_fig(gcf,[outdir 'Edgewise_Etasquared_' mat_names{m} '_' nameStart '.pdf']);
end
close('all');


end