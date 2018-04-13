function attribute_gifti_data_to_cifti(data_L, data_R, fname)
order_R = gifti(data_R);
order_L = gifti(data_L);
order_L = order_L.cdata;
order_R = order_R.cdata;

%medial_R = gifti('/data/cn4/haoxin/FP_Project/medial_wall.R.32k_fs_LR.func.gii');
medial_R = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');
%medial_L = gifti('/data/cn4/haoxin/FP_Project/medial_wall.L.32k_fs_LR.func.gii');
medial_L = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');
medial_R = medial_R.cdata;
medial_L = medial_L.cdata;


cifti_order = [order_L(~medial_L); order_R(~medial_R)];
cifti_order(size(cifti_order,1):66697,1) = 0;

cifti_write_wHDR(cifti_order, [], fname, 'dtseries')

