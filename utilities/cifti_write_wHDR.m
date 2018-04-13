function cifti_write_wHDR(data,ciftitemplatefile,outnamestem,varargin)  
%cifti_write_wHDR(data,[ciftitemplatefile],outnamestem,filetype)  


%where workbench is installed
workbenchdir = '/data/cn5/caterina/workbench/bin_linux64/';

if ~isempty(varargin)
    filetype = varargin{1};
else
    filetype = 'dtseries';
end
extension = ['.' filetype '.nii'];

if length(outnamestem) > length(extension) && strcmp(outnamestem(end-(length(extension)-1):end),extension)
    outnamestem = outnamestem(1:end-length(extension));
end

deletestuff= 0;
if isempty(ciftitemplatefile)
    ciftitemplatefile = '/data/cn4/evan/Scripts/cifti_template.func.gii';
elseif strcmp(ciftitemplatefile(end-(length(extension)-1):end),extension)
    system([workbenchdir 'wb_command -cifti-convert -to-gifti-ext ' ciftitemplatefile ' WritingTemp.func.gii'])
    ciftitemplatefile = 'WritingTemp.func.gii';
    deletestuff = 1;
end

dir = pwd;
save(gifti(single(data')),[outnamestem '.func.gii'],'ExternalFileBinary')

%Read template file and write out with correct header
bufsize = 524288;
ciftiheadertext = textread(ciftitemplatefile,'%s','delimiter','\r','bufsize',bufsize);
if deletestuff
    delete('WritingTemp.func*')
end
delete([outnamestem '.func.gii']);
fid = fopen([outnamestem '.func.gii'],'at'); %open the output file for writing
fclose(fid);


switch filetype
    
    case 'dconn'
        
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=70 && strcmp(thisline(1:70),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS"')
                thisline = '<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS" AppliesToMatrixDimension="0,1">';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=69 && strcmp(thisline(1:69),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_TIME_POINTS"')
            else
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            end           
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.dconn.nii'])
        
     case 'pconn'
        
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=70 && strcmp(thisline(1:70),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS"')
                thisline = '<MatrixIndicesMap AppliesToMatrixDimension="0,1" IndicesMapToDataType="CIFTI_INDEX_TYPE_PARCELS">';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=69 && strcmp(thisline(1:69),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_TIME_POINTS"')
            else
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            end           
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.pconn.nii'])

    case 'dtseries'
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
            end
            
            dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.dtseries.nii'])

end

delete([outnamestem '.func.*'])

