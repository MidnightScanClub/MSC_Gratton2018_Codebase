function [watertime water_corrmat water_ts_concat water_corrmat_concat tmask_all] = make_watershed_parcel_timecourse_cifti_func_TASK(tmasklist,timecoursedir,timestem,watershed_LR,outputdir,varargin)
% function [watertime water_corrmat water_ts_concat water_corrmat_concat tmask_all] = make_watershed_parcel_timecourse_cifti_func_TASK(tmasklist,timecoursedir,timestem,watershed_LR,outputdir,varargin)
% Watershed-based timecourse extraction, Left and Right
% TOL created
% CG edited to work for task data

% get information about dataset
[path name ext] = fileparts(watershed_LR);
cd(path)
if length(varargin) > 0
    [tempdataset tmasks subjects] = textread(tmasklist,'%s%s%s');
else
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

% read in CIFTI file
watershed = ft_read_cifti_mod(watershed_LR);
watershed = watershed.data(1:59412,:); % CG: remove medial wall info
waternum = unique(watershed); % get indices for each parcel
waternum(waternum==0) = [];
cd(outputdir)

% loop over subjects
for s = 1:length(subjects)
    
    
    disp(['Processing subject #' num2str(s) ': ' subjects{s}])
    tic
    timename = [subjects{s} '_' timestem];
    
    fname = [timecoursedir '/' timename '.dtseries.nii'];
    if ~exist(fname)
        watertime{s}(1,1:length(waternum)) = nan;
        goodsubs(s) = 0;
        disp(['missing data: ' timename]);
    else
        goodsubs(s) = 1;
        surf_timecourse = ft_read_cifti_mod(fname);
        surf_timecourse = surf_timecourse.data(1:59412,:);
        
        % loop over parcels
        for i = 1:length(waternum)
            
            % identify parcel indices and take the mean timeseries
            waterind = find(watershed==waternum(i));
            watertime{s}(:,i) = nanmean(surf_timecourse(waterind,:))'; %% added nanmean to address nan timeseries
            
            % show how many nan's there were to keep track of errors
            temp = surf_timecourse(waterind,:);
            if sum(isnan(temp(:))) > 0
                prop_nans = sum(isnan(temp(:)))/(length(temp(:)));
                disp(['Nans found for Sess=' num2str(s) ' in ROI=' num2str(i) ', ' num2str(prop_nans)]);
            end
            
            
        end
    end
    
end

% loop through subjects again
water_ts_concat = [];
for s = 1:length(subjects)
    if tmasks{s}(2) == 'n' %still has net at the start - deal with naming convention issue
        tmask_fname = ['/data' tmasks{s}(5:end)];
    else
        tmask_fname = tmasks{s};
    end
    
    % if no good data, fill with nans; otherwise, run correlations per session
    if goodsubs(s) == 0
        water_corrmat(:,:,s) = ones(length(waternum),length(waternum))*nan;
        tmask_all{s} = [0];
    else
        tmask = load(tmask_fname);
        water_corrmat(:,:,s) = FisherTransform(paircorr_mod(watertime{s}(logical(tmask),:)));
        water_ts_concat = [water_ts_concat; watertime{s}(logical(tmask),:)];
        tmask_all{s} = tmask;
    end
end

% also run correlations across concatenated data (in my experience, produces almost identical results)
water_corrmat_concat = FisherTransform(paircorr_mod(water_ts_concat));