function atlas_params = alter_atlas_params_removeUA(atlas_params_orig)
% assumes UA is the first row

good_ROIs = ones(1,atlas_params_orig.num_rois);
good_ROIs(atlas_params_orig.mods{1}) = 0;
good_ROIs = logical(good_ROIs);

atlas_params.num_rois = sum(good_ROIs);
atlas_params.dist_thresh = 20;

assignments = atlas_params_orig.mods_array(good_ROIs);
assignments_new = zeros(size(assignments));
IDs = unique(assignments);
for i = 1:length(IDs)
    mods{i} = find(assignments == IDs(i));
    assignments_new(assignments == IDs(i)) = i; %renumber so things come out correctly in later code
end
atlas_params.mods_array = assignments_new;
atlas_params.mods = mods;

[communities atlas_params.sorti] = sort(atlas_params.mods_array);
atlas_params.transitions = find(communities(1:end-1) ~= communities(2:end));
transitions_plusends = [1 atlas_params.transitions(:)' length(communities)];
atlas_params.centers = transitions_plusends(1:end-1) + ((transitions_plusends(2:end) - transitions_plusends(1:end-1))/2);

atlas_params.networks = atlas_params_orig.networks(2:end); %atlas_params_orig.networks(IDs); % a cludge, since ID #s don't match networks
atlas_params.colors = atlas_params_orig.colors(2:end,:); %atlas_params_orig.colors(IDs,:); %same as above
atlas_params.good_ROIs = good_ROIs;
%atlas_params.good_ROIs_init = good_ROIs_init;
atlas_params.notes = 'removed unlabeled ROIs';
atlas_params.atlas = [atlas_params_orig.atlas '_removedUA'];



end