function FramesPerCond_MSC(subject,task)
% function FramesPerBlock()
%
% Script to look at the amount of data contained per condition when
% separated into tasks, removing fixation periods
%
% subject = string with subject name (e.g., 'MSC01')
% task = string with task name (e.g., 'motor')
%
% CGratton - 6/3/2014
% CG - edited to work for MSC, 5/9/2016

% MAIN VARIABLES - TO CHANGE -
numTRs_skip = 4; % # of TRs to skip at the start of each run
snipsz = 2; % how short can segments be to be included (will actually be +1 from this number)

% dir info
topdir = '/data/nil-bluearc/GMT/Caterina/';
dsg_dir = [topdir 'SurfTask_analysis/' subject '/' task '_vol/design_matrices/']; % where design matrices are stored
FCDir = [topdir 'TaskFC/FCProc_' subject '_' task '_pass2/'];
FCDir_prescrub = [topdir 'TaskFC/FCProc_' subject '_' task '_pass1/'];

%array of dsg matrix conditions to keep
[dsg_conds_init,cond_types_init,cond_inds_init,cond_inds_noover_init] = get_design_condsANDinds(task);
if strcmp(subject,'MSC10') && strcmp(task,'mem')
    % to deal with session 6 not having scene conditions
    [dsg_conds_sess6,cond_types_sess6,cond_inds_sess6,cond_inds_noover_sess6] = get_design_condsANDinds_MSC10_mem(task);
end

% list of subjects:
%tmasklist = [FCDir_prescrub '/COHORTSELECT/NEW_CUT_TMASKLIST.txt'];
tmasklist = [FCDir_prescrub '/COHORTSELECT/NEW_TMASKLIST.txt']; % get all conditions, regardless if they have enough data
[sessions tmasks] = textread(tmasklist,'%s%s');

%load QC file
load([FCDir_prescrub '/QC.mat']); QC_prescrub = QC; clear QC;
load([FCDir '/QC.mat']);

%initialize file and write out a header
%fname = [FCDir '/frames_per_condition.txt'];
%fid = fopen(fname,'w+');
%fprintf(fid,['subnum\tsubID\tTOTframesOrig\tTOTframesFinal\t',...
%    'framesAllOrig\tnframesAllFinal\t',...
%    'framesTongOrig\tframesTongFinal\t'); %etc., etc., - make flexible by task


%loop through subjects
for s = 1:length(sessions)
    
    % select the right design matrix info (only an issue for MSC10, mem)
    if strcmp(subject,'MSC10') && strcmp(task,'mem') && strcmp(sessions{s},'vc39619A')
        dsg_conds = dsg_conds_sess6;
        cond_types = cond_types_sess6;
        cond_inds = cond_inds_sess6;
        cond_inds_noover = cond_inds_noover_sess6;
    else
        dsg_conds = dsg_conds_init;
        cond_types = cond_types_init;
        cond_inds = cond_inds_init;
        cond_inds_noover = cond_inds_noover_init;
    end
    
    % load the tmask
    tmask = load(tmasks{s});
    
    % load the design matrix
    fname = [dsg_dir subject '_' task '_' sessions{s} '_dsgn_matrix.txt'];
    dsg_data_orig = importdata(fname);
    
    % find the right QC index
    qc_id_prescrub = find_sub_indexQC(QC_prescrub,sessions{s});
    qc_id = find_sub_indexQC(QC,sessions{s});
    
    if qc_id
        %%% Deal with missing runs
        % get a list of the good runs and their run borders
        clear good_runs;
        for r = 1:length(QC(qc_id).restruns)
            good_runs(r) = str2num(QC(qc_id).restruns{r});
        end
        if strcmp(subject, 'MSC09') && strcmp(task,'motor') && strcmp(sessions{s},'vc39548') %EXCEPTION
            good_runs = good_runs - 1; %subtract 1, since runs are 2 and 3
        end
        
        %0 out data from missing blocks
        kept_run_mask = ones(1,QC_prescrub(qc_id_prescrub).numTRs);
        kept_run_mask = remove_run_mask(kept_run_mask,good_runs,QC_prescrub(qc_id_prescrub).runborders);
        if length(good_runs) == 1
            disp('one run');
        end
        
        %%%shorten this mask based on the first N timepoints in each run
        %since these aren't included in the dsgn matrix
        dsgmask_runstarts = create_run_mask(QC_prescrub(qc_id_prescrub).runborders,numTRs_skip,QC_prescrub(qc_id_prescrub).numTRs);
        kept_run_mask_dsgn = kept_run_mask(logical(dsgmask_runstarts));
        dsgmask_runstarts_short = create_run_mask(QC(qc_id).runborders,numTRs_skip,QC(qc_id).numTRs);
        
        %%% reformat design matrix
        dsgn_data_final = zeros(length(tmask),size(dsg_data_orig.data,2));
        dsgn_data_final_tmask = zeros(size(dsgn_data_final));
        dsgn_data_final_tmask_snip = zeros(size(dsgn_data_final));
        for col = 1:size(dsg_data_orig.data,2)
            temp = dsg_data_orig.data(logical(kept_run_mask_dsgn),col);
            dsgn_data_final(logical(dsgmask_runstarts_short),col) = temp;
            dsgn_data_final_tmask(:,col) = dsgn_data_final(:,col).*tmask;
            dsgn_data_final_tmask_snip(:,col) = sniptinymask(dsgn_data_final_tmask(:,col)>0.5,snipsz);
        end
        
        % make a plot to check the results
        if strcmp(task,'motor')
            figure;
            subplot(1,3,1)
            imagesc(dsgn_data_final>0.5,[-1 1]);
            title('conds');
            set(gca,'XTick',[1:length(dsg_conds)],'XTicklabel',dsg_conds);
            xticklabel_rotate();
            colormap('gray');
            subplot(1,3,2)
            imagesc(dsgn_data_final_tmask>0.5,[-1 1]);
            colormap('gray');
            title('conds, masked w tmask');
            set(gca,'XTick',[1:length(dsg_conds)],'XTicklabel',dsg_conds);
            xticklabel_rotate();
            subplot(1,3,3)
            imagesc(dsgn_data_final_tmask_snip,[-1 1]);
            colormap('gray');
            title('add contig th');
            set(gca,'XTick',[1:length(dsg_conds)],'XTicklabel',dsg_conds);
            xticklabel_rotate();
            fig_name = [FCDir subject '_' sessions{s} '_taskconds.png'];
            save_fig(gcf,fig_name);
        end
        
        % Find indices for each condition, w and w/o masking
        % Also count their lengths.
        for c = 1:length(cond_types)
            [TIndTot(s).(cond_types{c}), TIndFin(s).(cond_types{c})] = select_cond_indices(dsgn_data_final,dsgn_data_final_tmask,...
                cond_inds.(cond_types{c}),cond_inds_noover.(cond_types{c}),snipsz);
            
            %TIndTot(s).(cond_types{c}) = sum(dsgn_data_final(:,cond_inds.(cond_types{c}))>0.5,2)>0;
            %temp = sum(dsgn_data_final_tmask(:,cond_inds.(cond_types{c}))>0.5,2)>0;
            %TIndFin(s).(cond_types{c}) = sniptinymask(temp,snipsz); %redo here for final, since conditions grouped
            
            FramesTot(s).(cond_types{c}) = sum(TIndTot(s).(cond_types{c}));
            FramesFin(s).(cond_types{c}) = sum(TIndFin(s).(cond_types{c}));
            
            FDmeantot(s).(cond_types{c}) = mean(QC(qc_id).FD(logical(TIndTot(s).(cond_types{c}))));
            FDmeanfin(s).(cond_types{c}) = mean(QC(qc_id).FD(logical(TIndFin(s).(cond_types{c}))));
            FDtot(s).(cond_types{c}) = QC(qc_id).FD(logical(TIndTot(s).(cond_types{c})));
            FDfin(s).(cond_types{c}) = QC(qc_id).FD(logical(TIndFin(s).(cond_types{c})));
        end
    else
        % fill in empty values
        for c = 1:length(cond_types)
            TIndTot(s).(cond_types{c}) = 0;
            TIndFin(s).(cond_types{c}) = 0;
            FramesTot(s).(cond_types{c}) = sum(TIndTot(s).(cond_types{c}));
            FramesFin(s).(cond_types{c}) = sum(TIndFin(s).(cond_types{c}));
            FDmeantot(s).(cond_types{c}) = nan;
            FDmeanfin(s).(cond_types{c}) = nan;
            FDtot(s).(cond_types{c}) = nan;
            FDfin(s).(cond_types{c}) = nan;
        end
    end
    
end


%save out indices from each condition for use in fcimage_operations
fts_fname = [FCDir '/condindices.mat'];
save(fts_fname,'TIndTot','TIndFin','FramesTot','FramesFin','FDtot','FDfin','FDmeantot','FDmeanfin','dsg_conds','cond_types','cond_inds');


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = create_run_mask(run_borders,N,mask_size)
%function that creates a mask for first N TRs in each run

mask = squeeze(ones(1,mask_size));

for r = 1:size(run_borders,1)
    mask(run_borders(r,2):run_borders(r,2)+N-1) = 0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = remove_run_mask(mask,KEPTruns,run_borders)

%for r = 1:size(run_borders,1)
%    if length(find(KEPTruns == run_borders(r,1))) == 0
%
%        %if the run is NOT kept, zero it out
%        mask(run_borders(r,2):run_borders(r,3)) = 0;
%    end
%end
mask_new = zeros(size(mask));
for r = 1:length(KEPTruns)
    mask_new(run_borders(KEPTruns(r),2):run_borders(KEPTruns(r),3)) = 1;
end
mask = mask_new.*mask;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function full_ts_index = convert_to_ts_index(cond_ts,dsgts_mask)
% to convert back to full ts index, mainly just need to add back in TRs for
% each run (excluding removed runs)
% can do this by using original design mask

full_ts_index = zeros(size(dsgts_mask));
full_ts_index(dsgts_mask == 1) = cond_ts;
full_ts_index = logical(full_ts_index);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask2]=sniptinymask(mask1,snipsz)

tempmask1 = mask1;
[chunksize startnums endnums] = maskgaps(~tempmask1);
goodsize=endnums-startnums+1;
goodlost=goodsize<=snipsz;
for j=1:numel(goodlost)
    if goodlost(j)
        tempmask1(startnums(j):endnums(j))=0;
    end
end
mask2=tempmask1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chunksize startnums endnums] = maskgaps(mask)

% presumes that mask is 1=keep, 0=discard


if nnz(mask)==0 % if all discarded
    chunksize=-numel(mask);
elseif nnz(mask)==numel(mask) % if all kept
    chunksize=numel(mask);
else % if some kept and some discarded
    
    dmask=diff(mask);
    numchunks=nnz(dmask)+1;
    chunksize=ones(numchunks,1);
    chunknum=1;
    
    for i=1:numel(dmask)
        if dmask(i)==0
            chunksize(chunknum,1)=chunksize(chunknum,1)+1;
        elseif dmask(i)==1
            chunksize(chunknum,1)=-chunksize(chunknum,1);
            chunknum=chunknum+1;
        elseif dmask(i)==-1
            chunknum=chunknum+1;
        end
        
        if i==numel(dmask)
            if chunksize(chunknum-1,1)>0
                chunksize(chunknum,1)=-chunksize(chunknum,1);
            end
        end
        
    end
end

if ~isequal((sum(abs(chunksize))),numel(mask))
    disp('error');
end

endnums=cumsum(abs(chunksize));
endnums=endnums(chunksize<0);
startnums=endnums+chunksize(chunksize<0)+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dsg_conds,cond_types,cond_inds,cond_inds_noover] = get_design_condsANDinds(task)

switch task
    case 'motor'
        dsg_conds = {'Tongue','L_Hand','R_Hand','L_Leg','R_Leg','Trend','Trend','Baseline','Baseline'};
        cond_types = {'AllMotor','Tongue','L_Hand','R_Hand','L_Leg','R_Leg'};
        cond_inds.AllMotor = [1:5];
        cond_inds.Tongue = [1];
        cond_inds.L_Hand = [2];
        cond_inds.R_Hand = [3];
        cond_inds.L_Leg = [4];
        cond_inds.R_Leg = [5];
        
        % indices that shouldn't overlap with it
        cond_inds_noover.AllMotor = [];
        cond_inds_noover.Tongue = [2:5];
        cond_inds_noover.L_Hand = [1,3,4,5];
        cond_inds_noover.R_Hand = [1,2,4,5];
        cond_inds_noover.L_Leg = [1,2,3,5];
        cond_inds_noover.R_Leg = [1:4];
        
    case 'mixed'
        dsg_conds = {'Glass_startcue','Glass_sustained','Glass_endcue','Glass_coherent','Glass_random',...
            'NV_startcue','NV_sustained','NV_endcue','Noun','Verb','Trend','Trend','Baseline','Baseline'};
        %cond_types = {'AllGlass','Glass_Trials','AllSemantic','Semantic_Trials'};
        cond_types = {'AllGlass','AllSemantic'};
        cond_inds.AllGlass = [1:33];
        cond_inds.AllSemantic = [34:66];
        
        % indices that shouldn't overlap with it
        cond_inds_noover.AllGlass = cond_inds.AllSemantic;
        cond_inds_noover.AllSemantic = cond_inds.AllGlass;
        
        
    case 'mem'
        dsg_conds = {'Indoor_1st','Indoor_2nd','Indoor_3rd','Outdoor_1st','Outdoor_2nd','Outdoor_3rd',...
            'Female_1st','Female_2nd','Female_3rd','Male_1st','Male_2nd','Male_3rd',...
            'Abstract_1st','Abstract_2nd','Abstract_3rd','Concrete_1st','Concrete_2nd','Concrete_3rd',...
            'Trend','Trend','Trend','Baseline','Baseline','Baseline'};
        cond_types = {'AllMem','Scene','Face','Word','pres1','pres2','pres3'};
        cond_inds.AllMem = [1:144];
        cond_inds.Scene = [1:48];
        cond_inds.Face = [49:96];
        cond_inds.Word = [97:144];
        cond_inds.pres1 = [1:8,25:32,49:56,73:80,97:104,121:128];
        cond_inds.pres2 = [9:16,33:40,57:64,81:88,105:112,129:136];
        cond_inds.pres3 = [17:24,41:48,65:72,89:96,113:120,137:144];
        
        % indices that shouldn't overlap with it
        cond_inds_noover.AllMem = [];
        cond_inds_noover.Scene = [cond_inds.Face cond_inds.Word];
        cond_inds_noover.Face = [cond_inds.Scene cond_inds.Word];
        cond_inds_noover.Word = [cond_inds.Scene cond_inds.Face];
        cond_inds_noover.pres1 = [cond_inds.pres2 cond_inds.pres3];
        cond_inds_noover.pres2 = [cond_inds.pres1 cond_inds.pres3];
        cond_inds_noover.pres3 = [cond_inds.pres1 cond_inds.pres2];
        
end
end

function [dsg_conds,cond_types,cond_inds,cond_inds_noover] = get_design_condsANDinds_MSC10_mem(task)
%MSC10 sess 6 has no scenes
switch task
    case 'mem'
        dsg_conds = {'Female_1st','Female_2nd','Female_3rd','Male_1st','Male_2nd','Male_3rd',...
            'Abstract_1st','Abstract_2nd','Abstract_3rd','Concrete_1st','Concrete_2nd','Concrete_3rd',...
            'Trend','Trend','Trend','Baseline','Baseline','Baseline'};
        cond_types = {'AllMem','Scene','Face','Word','pres1','pres2','pres3'};
        cond_inds.AllMem = [1:96]; %[1:144];
        cond_inds.Scene = [];
        cond_inds.Face = [1:48]; %[49:96];
        cond_inds.Word = [49:96]; %[97:144];
        cond_inds.pres1 = [1:8,25:32,49:56,73:80]; %[1:8,25:32,49:56,73:80,97:104,121:128];
        cond_inds.pres2 = [9:16,33:40,57:64,81:88]; %[9:16,33:40,57:64,81:88,105:112,129:136];
        cond_inds.pres3 = [17:24,41:48,65:72,89:96]; %[17:24,41:48,65:72,89:96,113:120,137:144];
        
        % indices that shouldn't overlap with it
        cond_inds_noover.AllMem = [];
        cond_inds_noover.Scene = [cond_inds.Face cond_inds.Word];
        cond_inds_noover.Face = [cond_inds.Scene cond_inds.Word];
        cond_inds_noover.Word = [cond_inds.Scene cond_inds.Face];
        cond_inds_noover.pres1 = [cond_inds.pres2 cond_inds.pres3];
        cond_inds_noover.pres2 = [cond_inds.pres1 cond_inds.pres3];
        cond_inds_noover.pres3 = [cond_inds.pres1 cond_inds.pres2];
        
end
end

function [TIndTot,TIndFin] = select_cond_indices(dsgmat,dsgmat_tmask,inds,inds_noover,snipsz)

temp1 = sum(dsgmat(:,inds)>0.5,2)>0;
temp2 = sum(dsgmat(:,inds_noover)>0.5,2)>0;
%temp1 = sum(dsgmat(:,inds),2)>0;
%temp2 = sum(dsgmat(:,inds_noover),2)>0;
TIndTot = temp1 .* ~temp2;

temp1 = sum(dsgmat_tmask(:,inds)>0.5,2)>0;
temp3 = temp1 .* ~temp2;
TIndFin = sniptinymask(temp3,snipsz); %redo here for final, since conditions grouped

end