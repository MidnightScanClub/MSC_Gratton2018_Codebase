function h = figure_corrmat_network_generic(matrix_orig,atlas_params,varargin)
% function h = figure_corrmat_powernetwork(matrix_orig,atlas_params,varargin)
% Makes a correlation matrix, with lines marking the divisions between
% networks, and network labels on sides
%
% Originally from Tim - works for his network structure
% C Gratton - modified to make this work with any network structure, as
% specified in the atlas_params input
%
% Inputs:
%   matrix_orig: the matrix that will be plotted
%   atlas_params: a structure with information about the atlas (see atlas_parameters.m)
%   varargin: low,hi color limits

% information about networks and colors
networks = atlas_params.networks;
colors_new = atlas_params.colors;

% sort matrix into correct order
matrix = matrix_orig(atlas_params.sorti,atlas_params.sorti);

% open figure
h = figure('Color',[0.8275 0.8275 0.8275],'Position',[56 143 1295 807]); %[56 143 1095 807]

% plot out matrix
if nargin>2
    if ~isempty(varargin{1}) & ~isempty(varargin{2})
        climlow = varargin{1};
        climhigh = varargin{2};
        imagesc(matrix,[climlow climhigh]);
    else
        imagesc(matrix);
    end
else
    imagesc(matrix);
end

% My favorite colormap - edit if you want a different one
load /data/cn5/caterina/PDgrant/scripts/better_jet_colormap.mat;
colormap(better_jet_colormap_diff);

% put lines between the networks
vline_new(atlas_params.transitions,'k',3);
hline_new(atlas_params.transitions,'k',3);
tickpos = atlas_params.centers;
ax = axis;

% put ticks in the right places
set(gca,'XTick',tickpos,'Xlim',[ax(1) ax(2)]);

% take out the current tick labels
set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

% label the networks on the two axes
tx= text(tickpos,ones(1,length(tickpos))*(atlas_params.num_rois+1),networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

for i = 1:length(tx)
     set(tx(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
end
set(gca,'FontWeight','bold','FontSize',10);


ty= text(-1*ones(1,length(tickpos)),tickpos-5,networks);
set(ty,'HorizontalAlignment','right','VerticalAlignment','top')

for i = 1:length(ty)
    set(ty(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
end
colorbar;
set(gca,'FontWeight','bold','FontSize',10);

% if provided, give the plot a title
if nargin>4
    if ~isempty(varargin{4})
        title(titletext,'FontWeight','bold','FontSize',10);
    end
end

% make square
axis square;