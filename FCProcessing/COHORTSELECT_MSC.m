function COHORTSELECT(QCmat,datafile,vcidfile,varargin)

% jdp 1/13/13

outdir='COHORTSELECT';
if ~exist(outdir,'dir')
	mkdir(outdir);
end

% all the stuff from subjects
load(QCmat);
subset=1:numel(QC);

% load in the original datafile because I'm a dumbass
[df.filepath df.vcnum df.prmfile df.TR df.TRskip] = textread(datafile,'%s%s%s%f%d');

% read in the vcid file, incorporate into QC structure
for i=1:size(QC,2)
    dfile.prepstem{i,1}=QC(i).vcnum;
end

[vcid.vc vcid.age vcid.sex vcid.uniqueID vcid.twinstogetheruniqueID] = vcIDfilereader(dfile.prepstem,vcidfile);
for i=1:size(QC,2)
    QC(i).age=vcid.age(i,1);
    QC(i).sex=vcid.sex{i,1};
    QC(i).uniqueID=vcid.uniqueID(i,1);
    QC(i).twinstogetheruniqueID=vcid.twinstogetheruniqueID(i,1);
end

% get temporal masks if existing
if ~isempty(varargin)
	[vcs tmaskfile] = textread(varargin{1,1},'%s%s');
	if ~isequal(vcs,dfile.prepstem)
		error('VC ordering in datafile and tmasklist are not the same.');
	end
end

processes=QC(1).process';
stages=numel(processes);

% initialize masks
threshold_FD=Inf;
threshold_DV=repmat(Inf,[1 stages]);
threshold_SD=repmat(Inf,[1 stages]);
augment_FD=[0; 0];
augment_DV=repmat([0; 0],[1 stages]);
augment_SD=repmat([0; 0],[1 stages]);

% set thresholds
threshold_FD=input('What FD threshold? (0.2 recommended) :');
augment_FD(1)=input('Augment FD how far forward (0 frames recommended with GSR) :');
augment_FD(2)=input('Augment FD how far back (0 frames recommended with GSR) :');
snipsz=input('What minimum number of contiguous TRs needed? (5 recommended) :');
runmin=input('What minimum number of volumes in run needed? (50 recommended) :');
totalmin=input('What minimum number of volumes from all runs is needed? (100 recommended) :');
% threshold_DV(2)=20;

save([outdir '/cohortselect_MSC_params.mat' ],...
    'threshold_FD','augment_FD','snipsz','runmin','totalmin');


%%%
%%% HISTOGRAMS OF ENTIRE GROUP
%%%
% 
% % image histograms
% [tot]=QCdistrib(QC,subset);
% 
% clf;
% [fdhist fdlimz]=hist(tot.FD,2000);
% h1=bar(fdlimz,fdhist);
% set(h1,'facecolor','r','edgecolor','r');
% hold('on');
% z=cumsum(fdhist);
% ylimz=get(gca,'ylim');
% z=z/max(z)*ylimz(2);
% plot(fdlimz,z,'k');
% hold('off');
% xlim([0 1]);
% vline(threshold_FD,'k');
% saveas(gcf,[outdir '/FDhist.tiff'],'tiff');
% 
% for k=1:stages
%     if ~isinf(threshold_DV(1,k))
%         clf;
%         [fdhist fdlimz]=hist(tot.DV(:,k),2000);
%         h1=bar(fdlimz,fdhist);
%         set(h1,'facecolor','b','edgecolor','b');
%         hold('on');
%         z=cumsum(fdhist);
%         ylimz=get(gca,'ylim');
%         z=z/max(z)*ylimz(2);
%         plot(fdlimz,z,'k');
%         hold('off');
%         xlim([0 25]);
%         vline(threshold_DV(1,k),'k');
%         xlim([tot.DV_median(1,k)-2*tot.DV_std(1,k) tot.DV_median(1,k)+4*tot.DV_std(1,k)]);
% %         vline([tot.DV_median(1,k) tot.DV_median(1,k)-tot.DV_std(1,k) tot.DV_median(1,k)+tot.DV_std(1,k)],'w');
%         saveas(gcf,[outdir '/DVhist_QCstage' num2str(k) '.tiff'],'tiff');
%     end
% end
% 
% for k=1:stages
%     if ~isinf(threshold_SD(1,k))
%         clf;
%         [fdhist fdlimz]=hist(tot.SD(:,k),2000);
%         h1=bar(fdlimz,fdhist);
%         set(h1,'facecolor','g','edgecolor','g');
%         hold('on');
%         z=cumsum(fdhist);
%         ylimz=get(gca,'ylim');
%         z=z/max(z)*ylimz(2);
%         plot(fdlimz,z,'k');
%         hold('off');
%         xlim([0 25]);
%         vline(threshold_SD(1,k),'k');
%         xlim([tot.SD_median(1,k)-2*tot.SD_std(1,k) tot.SD_median(1,k)+4*tot.SD_std(1,k)]);
%         vline([tot.SD_median(1,k) tot.SD_median(1,k)-tot.SD_std(1,k) tot.SD_median(1,k)+tot.SD_std(1,k)],'w');
%         saveas(gcf,[ outdir '/SDhist_QCstage' num2str(k) '.tiff'],'tiff');
%     end
% end


for i=subset
    QC(i).MSK_skipmask=ones(size(QC(i).FD,1),1);
    for j=1:size(QC(i).runborders,1)
        QC(i).MSK_skipmask(QC(i).runborders(j,2):QC(i).runborders(j,2)+QC(i).TRskip-1)=0;
    end
    
    counter=0;
    if ~isinf(threshold_FD)
        counter=counter+1;
        QC(i).MSK_QC(:,counter)=expandmask(~(QC(i).FD>threshold_FD),augment_FD(1),augment_FD(2),QC(i).runborders);
    end
    
    for k=1:stages
        if ~isinf(threshold_DV(1,k))
            counter=counter+1;
            QC(i).MSK_QC(:,counter)=expandmask(~(QC(i).DV_GM>threshold_DV(k)),augment_DV(1,k),augment_DV(2,k),QC(i).runborders);
        end
    end

    for k=1:stages
        if ~isinf(threshold_SD(1,k))
            counter=counter+1;
            QC(i).MSK_QC(:,counter)=expandmask(~(QC(i).SD_GM>threshold_SD(k)),augment_SD(1,k),augment_SD(2,k),QC(i).runborders);
        end
    end
    
    % combine skipmask and QC masks (frame excluded if ANY indicates bad)
    QC(i).MSK_SKIP_TOTQC=QC(i).MSK_skipmask & ~sum(~(QC(i).MSK_QC),2);
    
    % trim out tiny bits (<X contiguous TRs)
    QC(i).MSK_SKIP_TOTQC_SNIP=sniptinymask(QC(i).MSK_SKIP_TOTQC,QC(i).runborders,snipsz);
    
    % zero out runs with <X volumes remaining
    [QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN borderqual QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN_new]=zerorun(QC(i).MSK_SKIP_TOTQC_SNIP,QC(i).runborders,runmin);
    QC(i).restruns_new=QC(i).restruns(~~borderqual);
    if sum(borderqual)~=numel(borderqual)
        fprintf('%d\t%s\n',i,QC(i).vcnum);
        % CG - also calculate and print number of volumes per run
        for r = 1:size(QC(i).runborders,1)
            good_frames_run(r) = sum(QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN(QC(i).runborders(r,2):QC(i).runborders(r,3)));
        end
        [borderqual QC(i).runborders good_frames_run']
    end
    
    remaining(i)=sum(QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN);
    disp(sprintf('Sub: %d\t%d\t%d\t%.03f',i,length(QC(i).MSK_QC),remaining(i),remaining(i)/length(QC(i).MSK_QC)));
    
end
    subsuse = remaining>totalmin;
subsuse_ind = find(subsuse);
disp(['Number of subjects remaining = ' num2str(sum(subsuse))])
for i=subset
     % get absolute displacement
    QC(i).sumtrans=sum(abs(QC(i).DTMVM(:,1:3)),2);
    QC(i).sumrot=sum(abs(QC(i).DTMVM(:,4:6)),2);
    
    clf;
    subplot(3,1,1);
    plot(QC(i).FD,'r','linewidth',2);
    hold on;
    plot(QC(i).sumtrans,'-.','color',[128 0 0]/255,'linewidth',1);
    plot(QC(i).sumrot,':','color',[128 0 0]/255,'linewidth',1);
    xlim([0 numel(QC(i).sumrot)]); ylim([0 1]); set(gca,'xticklabel','');
    vline(QC(i).runborders(:,3),'k');
    hline(threshold_FD,'k');
    
    subplot(3,1,2);
    plot(QC(i).DV_GM(:,2),'b','linewidth',2);
    xlim([0 numel(QC(i).sumrot)]); ylim([0 30]); set(gca,'xticklabel','');
    vline(QC(i).runborders(:,3),'k');
    hline(threshold_DV(1,2),'k');
    
    subplot(3,1,3);
    imagesc([QC(i).MSK_skipmask QC(i).MSK_SKIP_TOTQC QC(i).MSK_SKIP_TOTQC_SNIP QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN]');
    saveas(gcf,[ outdir '/' num2str(i) '_' QC(i).vcnum '.tiff'],'tiff');
end



% write new prmfiles, and tmask files
for i=subset
    
    % write NEW prmfile
    QC(i).newprmfile=[ pwd '/' outdir '/' QC(i).vcnum '_NEW.prm' ];
    QC(i).newtmaskfile=[ pwd '/' outdir '/' QC(i).vcnum '_NEW_TMASK.txt' ];
    QC(i).oldtmaskfile=[ pwd '/' outdir '/' QC(i).vcnum '_OLDCOMPATIBLE_TMASK.txt' ];
    
    fid=fopen(QC(i).newprmfile,'w');
    fprintf(fid,'set boldruns = (');
    for j=1:numel(QC(i).restruns_new)
        %fprintf(fid,'%d ',cell2mat(QC(i).restruns_new(j))); % CG - changed here; didn't seem to work properly before?
        fprintf(fid,'%s ',cell2mat(QC(i).restruns_new(j)));
    end
    fprintf(fid,')');
    fclose(fid);
    
    dlmwrite(QC(i).newtmaskfile,QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN_new,'\t');
    dlmwrite(QC(i).oldtmaskfile,QC(i).MSK_SKIP_TOTQC_SNIP_RUNMIN,'\t');

end

% new tmasklists and datalist
fid1=fopen([ pwd '/' outdir '/NEW_TMASKLIST.txt' ],'w');
fid2=fopen([ pwd '/' outdir '/OLD_TMASKLIST.txt' ],'w');
fid3=fopen([ pwd '/' outdir '/NEW_DATALIST.txt' ],'w');
for i=subset
    fprintf(fid1,'%s\t%s\n',QC(i).vcnum,QC(i).newtmaskfile);
    fprintf(fid2,'%s\t%s\n',QC(i).vcnum,QC(i).oldtmaskfile);
    fprintf(fid3,'%s\t%s\t%s\t%f\t%d\n',df.filepath{i,1},QC(i).vcnum,QC(i).newprmfile,QC(i).TR,QC(i).TRskip);
end
fclose(fid1);
fclose(fid2);
fclose(fid3);

% Write out datalist and tmasklist with subjects cut
fid1=fopen([ pwd '/' outdir '/NEW_CUT_TMASKLIST.txt' ],'w');
fid2=fopen([ pwd '/' outdir '/NEW_CUT_DATALIST.txt' ],'w');
for i=1:length(subsuse_ind)
    ind = subsuse_ind(i);
    fprintf(fid1,'%s\t%s\n',QC(ind).vcnum,QC(ind).newtmaskfile);
    fprintf(fid2,'%s\t%s\t%s\t%f\t%d\n',df.filepath{ind,1},QC(ind).vcnum,QC(ind).newprmfile,QC(ind).TR,QC(ind).TRskip);
end
fclose(fid1);
fclose(fid2);

disp hi;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask1 borderqual mask2]=zerorun(mask1,borders,runmin)

mask1=double(mask1);
r=size(borders,1);
for j=1:r
    borderqual(j,1)=1;
    if nnz(mask1(borders(j,2):borders(j,3)))<runmin
        borderqual(j,1)=0;
        mask1(borders(j,2):borders(j,3),1)=-1;
    end
end

mask2=mask1;
mask2(mask2==-1)=[];
mask1(mask1==-1)=0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask2]=sniptinymask(mask1,borders,snipsz)

runs=size(borders,1);
for r=1:runs
    tempmask1=mask1(borders(r,2):borders(r,3));
    [chunksize startnums endnums] = maskgaps(~tempmask1);
    goodsize=endnums-startnums+1;
    goodlost=goodsize<=snipsz;
    for j=1:numel(goodlost)
        if goodlost(j)
            tempmask1(startnums(j):endnums(j))=0;
        end
    end
    mask2(borders(r,2):borders(r,3))=tempmask1;
    
end
    
mask2=mask2';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask2]=expandmask(mask1,bnum,fnum,borders)

runs=size(borders,1);
for r=1:runs
    tempmask1=mask1(borders(r,2):borders(r,3));
    tempmask2=mask1(borders(r,2):borders(r,3));
    for x=1:numel(tempmask1)
        if tempmask1(x)==0
            for j=1:bnum
                if x-j>0
                    tempmask2(x-j)=0;
                end
            end
            for j=1:fnum
                if x+j<=numel(tempmask1)
                    tempmask2(x+j)=0;
                end
            end
        end
    end
    
    mask2(borders(r,2):borders(r,3))=tempmask2;
end
mask2=mask2';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tot]=QCdistrib(QC,subset)

stages=size(QC(1).DV_GM,2);

temptot=[];
for i=subset
    tt=~~QC(i).filtertmask;
    temptot=[temptot; QC(i).FD(tt)];
end
tot.FD=temptot;
tot.FD_mean=mean(tot.FD);
tot.FD_median=median(tot.FD);
tot.FD_std=std(tot.FD);
tot.FD_sort=sort(tot.FD);
tot.FD_cumsum=cumsum(tot.FD_sort)/sum(tot.FD_sort);

for j=1:stages
    temptot=[];
    for i=subset
        tt=~~QC(i).filtertmask;
        temptot=[temptot; QC(i).DV_GM(tt,j)];
    end
    tot.DV(:,j)=temptot;
    tot.DV_mean(1,j)=mean(tot.DV(:,j));
    tot.DV_median(1,j)=median(tot.DV(:,j));
    tot.DV_std(1,j)=std(tot.DV(:,j));
    tot.DV_sort(:,j)=sort(tot.DV(:,j));
    tot.DV_cumsum(:,j)=cumsum(tot.DV_sort(:,j))/sum(tot.DV_sort(:,j));
end

for j=1:stages
    temptot=[];
    for i=subset
        tt=~~QC(i).filtertmask;
        temptot=[temptot; QC(i).SD_GM(tt,j)];
    end
    tot.SD(:,j)=temptot;
    tot.SD_mean(1,j)=mean(tot.SD(:,j));
    tot.SD_median(1,j)=median(tot.SD(:,j));
    tot.SD_std(1,j)=std(tot.SD(:,j));
    tot.SD_sort(:,j)=sort(tot.SD(:,j));
    tot.SD_cumsum(:,j)=cumsum(tot.SD_sort(:,j))/sum(tot.SD_sort(:,j));
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


