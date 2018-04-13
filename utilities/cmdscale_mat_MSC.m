function [Y E] = cmdscale_mat_MSC(mat,groups,dist_type,colors,varargin)
% This function performs multi-dimensional scaling on input matrices and
% displays a 2-dimensional plot of the relative positions of the input
% matrices colored by group. MDS is performed using euclidean distance.
% Inputs:
% mat - node x node x subject array of matrices from all subjects
% groups - identifies which subject belongs to which group
% dist_type - type of distance calculation (e.g., euclidan or hamming)
% varargin -
%   (1) group names can be specified by a third input, e.g.
%       {'group1';'group2'} and will be displayed in a legend.
% Outputs:
% Y - configuration matrix
% E - eigenvalues of Y*Y'
%
% TOL, 09/14
% CG edits:
%   took out figure plotting command; do this outside this script
%   changed group indexing so it can deal with nonconsecutively numbered groups
% CG2: edited to do input dimensions

numnodes = size(mat,1);
group_vals = unique(groups); % CG - added this so we can do non consecutive groups
numgroups = length(group_vals); 
%colors = distinguishable_colors(numgroups);
if nargin > 4
    if length(varargin{1}) > 0
        dolabels = 1;
        groupnames = varargin{1};
    else
        dolabels = 0;
    end
end
if nargin > 5
    dims = varargin{2};
else
    dims = [1,2,3];
end

mask = ones(numnodes);
mask = triu(mask,1);
for s = 1:size(mat,3)
    temp = mat(:,:,s);
    mat_col(:,s) = temp(logical(mask));
end

% Multi-dimensional scaling
D = pdist(double(mat_col'),dist_type);
[Y E] = cmdscale(D);

% Display result
%figure('Color','white')    % CG - assumes figure has already been made outside
hold;
for c = 1:numgroups
    inds = groups==group_vals(c); % CG - changed this to do non consecutve groups
    %plot(Y(inds,1),Y(inds,2),'.','Color',colors(c,:),'MarkerSize',20);
        plot3(Y(inds,dims(1)),Y(inds,dims(2)),Y(inds,dims(3)),'.','Color',colors(c,:),'MarkerSize',20)
end
xlabel(['dim' num2str(dims(1))]); 
ylabel(['dim' num2str(dims(2))]); 
zlabel(['dim' num2str(dims(3))]); 

set(gca,'FontWeight','bold','FontSize',14);
if dolabels
    legend(groupnames,'FontWeight','bold','FontSize',14);
end