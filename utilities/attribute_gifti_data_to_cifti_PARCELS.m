function attribute_gifti_data_to_cifti_PARCELS(data_L, data_R, fname)
order_R = gifti(data_R);
order_L = gifti(data_L);
order_L = order_L.cdata;
order_R = order_R.cdata; %+ 33000; %so no chance that the parcel numbers will overlap
%order_R(order_R==33000) = 0; %but put 0 at 0

% renumber parcels within each
final_order_R = zeros(size(order_R));
parcel_R_list = unique(order_R);
parcel_R_list(parcel_R_list == 0) = []; % remove 0
for p = 1:length(parcel_R_list)
    final_order_R(order_R == parcel_R_list(p)) = p;
end
final_order_L = zeros(size(order_L));
parcel_L_list = unique(order_L);
parcel_L_list(parcel_L_list == 0) = []; % remove 0
for p = 1:length(parcel_L_list)
    final_order_L(order_L == parcel_L_list(p)) = p+length(parcel_R_list);
end

%medial_R = gifti('/data/cn4/haoxin/FP_Project/medial_wall.R.32k_fs_LR.func.gii');
medial_R = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');
%medial_L = gifti('/data/cn4/haoxin/FP_Project/medial_wall.L.32k_fs_LR.func.gii');
medial_L = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');
medial_R = medial_R.cdata;
medial_L = medial_L.cdata;

cifti_order = [final_order_L(~medial_L); final_order_R(~medial_R)];
cifti_order(size(cifti_order,1):66697,1) = 0;

cifti_write_wHDR(cifti_order, [], fname, 'dtseries')

