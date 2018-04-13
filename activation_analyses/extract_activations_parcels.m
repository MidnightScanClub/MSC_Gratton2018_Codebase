function [watertime water_ts_concat] = extract_activations_parcels(watershed_LR,vcid_list,actDir,startname,endname)
% function [watertime water_ts_concat] = extract_activations_parcels(watershed_LR,vcid_list,actDir,startname,endname)
% Watershed-based timecourse extraction, Left and Right
% CG - altered to work on activaiton maps

[path name ext] = fileparts(watershed_LR);

watershed = ft_read_cifti_mod(watershed_LR);
%watershed = watershed.data;
watershed = watershed.data(1:59412,:); % take just cortex
waternum = unique(watershed);
waternum(waternum==0) = [];

for s = 1:length(vcid_list)
    
    
    disp(['Processing subject #' num2str(s) ': ' vcid_list{s}])
    
    if strcmp(vcid_list{s},'vc39619A') && strcmp(endname,'AllScene')
        watertime{s} = ones(size(watertime{s-1}))*nan;
    else
        
        surf_timecourse = ft_read_cifti_mod([actDir '/' startname '_' vcid_list{s} '_' endname '.dscalar.nii']);
        surf_timecourse = surf_timecourse.data(1:59412,:); % take just cortex
        for i = 1:length(waternum)
            waterind = find(watershed==waternum(i));
            watertime{s}(:,i) = nanmean(surf_timecourse(waterind,:))'; 
            % and show how many nan's there were
            temp = surf_timecourse(waterind,:);
            if sum(isnan(temp(:))) > 0
                prop_nans = sum(isnan(temp(:)))/(length(temp(:)));
                disp(['Nans found for Sess=' num2str(s) ' in ROI=' num2str(i) ', ' num2str(prop_nans)]);
            end
        end
    end

end

water_ts_concat = [];
for s = 1:length(vcid_list) 
    water_ts_concat = [water_ts_concat; watertime{s}];
end
