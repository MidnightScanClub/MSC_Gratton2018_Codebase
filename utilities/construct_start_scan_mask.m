function start_scan_mask = construct_start_scan_mask(QC,start_scan_cut)

for i = 1:length(QC)
    start_scan_mask{i} = ones(QC(i).runborders(end,end),1);
    for b = 1:size(QC(i).runborders,1)
        start_scan_mask{i}(QC(i).runborders(b,2):QC(i).runborders(b,2)+start_scan_cut-1) = 0;
    end
    
end
end