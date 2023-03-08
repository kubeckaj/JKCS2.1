function [fluxes,clust_bound] = track_fluxes(Ca,Cb,basename,varargin)
% Find the main flux chains through the system for a 2-component case
% Input:
% Ca, Cb = acid and base concentrations (cm^-3 or ppt, see C_mon_notation below)
% basename = name of the base molecule, e.g. 'D'
% + optional input (see below)
% Output:
% fluxes = a structure containing the fluxes (cm^-3 s^-1), and the names of the starting and ending clusters
% clust_bound = clusters from which the outgoing fluxes originate

addpath ./functions -end


%%%%%%%%% Some settings (change if needed) %%%%%%%%%

% Charge: 1=neutral, 2=neg., 3= pos.
ch=1; % Neutral is the default

acidname='A';
acidname_neg='B';

name_mon={acidname,basename};
C_mon=[Ca,Cb];
C_mon_notation={'e','f'};

suffix='';

% Default criterion for plotting the flux to a certain cluster or out
% (fraction of the total flux into the cluster/out)
crit_out=0.05;
crit_clust=0.05;

% See what is given as the function input
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'ch')
            ch=varargin{i+1};
        elseif strcmpi(varargin{i},'C_mon_notation') | strcmpi(varargin{i},'notation')
            C_mon_notation=varargin{i+1};
        elseif strcmpi(varargin{i},'suffix')
            suffix=varargin{i+1};
        elseif strcmpi(varargin{i},'dirpath') | strcmpi(varargin{i},'dp')
            dirpath=varargin{i+1};
        elseif strcmpi(varargin{i},'crit_out')
            crit_out=varargin{i+1};
        elseif strcmpi(varargin{i},'crit_clust')
            crit_clust=varargin{i+1};
        elseif strcmpi(varargin{i},'temp')
            temp=varargin{i+1};
        elseif strcmpi(varargin{i},'IPR')
            IPR=varargin{i+1};
        end
    end
end

% Default directory path, if one is not given
if ~exist('dirpath','var') && exist('temp','var')
    dirpath=['.\',cell2mat(name_mon),'_',num2str(temp),'\'];
end

%if strcmpi(basename,'N')
%    maxacids=5;
%    maxbases=5;
%elseif strcmpi(basename,'M')
%    maxacids=3;
%    maxbases=3;
%elseif strcmpi(basename,'D')
%    maxacids=4;
%    maxbases=4;
%elseif strcmpi(basename,'T')
%    maxacids=2;
%    maxbases=2;
%else
%    maxacids=4;
%    maxbases=4;
%end
maxacids=20
maxbases=0

% Number of colors in the arrow plot
nlevels=5;


%%%%%%%%% End of the settings section %%%%%%%%%

fluxes=struct;


% First see what's coming out
if exist('IPR','var')
    [~,~,roundoutname,roundout] = plotflux_out(name_mon,C_mon,'C_mon_notation',C_mon_notation,...
        'crit',crit_out,'dp',dirpath,'suffix',suffix,'IPR',IPR);
else
    [~,~,roundoutname,roundout] = plotflux_out(name_mon,C_mon,'C_mon_notation',C_mon_notation,...
        'crit',crit_out,'dp',dirpath,'suffix',suffix);
end
    
k=0;
for i=1:length(roundout)
    lvalid=0;
    if ~isempty(strfind(roundoutname{i},'+'))
        clusters=regexp(roundoutname{i},'\+','split');        
        for j=1:2
            clusters_ch(j)=solve_charge(clusters{j});
        end
        % Include this flux if one of the colliders is of the wanted charge,
        % unless neutral pathways are wanted and the neutral collider is the smaller one
        if ismember(ch,clusters_ch)
            if ch==1
                [largest,smallest,maxmols,minmols]=max_mols(clusters{1},clusters{2});
                nl=strcmp(clusters,largest);
                ns=strcmp(clusters,smallest);
                if clusters_ch(nl)==ch
                    lvalid=1;
                    validclust=largest;
                elseif maxmols==minmols && clusters_ch(ns)==ch
                    lvalid=1;
                    validclust=smallest;
                end
            else
                lvalid=1;
                nc=clusters_ch==ch;
                validclust=clusters{nc};
           end
        end
        
        if lvalid==1
            combined=combine_clusters(clusters{1},clusters{2});
            k=k+1;
            fluxes(k).start=validclust;
            fluxes(k).end=combined;
            fluxes(k).val=roundout(i);
        end
    end
end

if k==0
    disp(['No contribution from ch = ',num2str(ch),' to the flux out for crit_out = ',num2str(crit_out)])
    return
end

% Total rate out of the system
tot_out=sum(roundout,'omitnan');
% Total rate out via the wanted charge
tot_out_ch=sum([fluxes.val],'omitnan');
titlestr=['Fraction of the total formation rate ',num2str(tot_out_ch/tot_out)];
% Another (optional) criterion for plotting the flux: threshold (compared to
% the total flux out along the major routes of the wanted charge)
%thflux=0.1*tot_out_ch;
thflux=0; % Set this to zero if you don't want to use it

% Save the names of the boundary clusters
clust_bound={fluxes.start};

nflux_out=k;
tracked={};
track_clus={};
for i1=1:nflux_out
    % Solve the origin of this, if it's a cluster
    if check_cluster(fluxes(i1).start) && ~ismember(fluxes(i1).start,tracked)
        start_ch=solve_charge(fluxes(i1).start);
        if start_ch==ch
            if exist('IPR','var')
                [~,~,mainfluxnames,mainfluxes,~]=plotflux(name_mon,C_mon,fluxes(i1).start,'to',0,...
                    'C_mon_notation',C_mon_notation,'crit',crit_clust,'dp',dirpath,'suffix',suffix,...
                    'lsame_ch','IPR',IPR);
            else
                [~,~,mainfluxnames,mainfluxes,~]=plotflux(name_mon,C_mon,fluxes(i1).start,'to',0,...
                    'C_mon_notation',C_mon_notation,'crit',crit_clust,'dp',dirpath,'suffix',suffix,...
                    'lsame_ch'); 
            end
            for i2=1:length(mainfluxes{1})
                if mainfluxes{1}(i2)>=thflux && ~strcmpi(mainfluxnames{1}{i2},'others')
                    k=k+1;
                    fluxes(k).end=fluxes(i1).start;
                    fluxes(k).start=mainfluxnames{1}{i2};
                    fluxes(k).val=mainfluxes{1}(i2);
                    if ~ismember(mainfluxnames{1}{i2},track_clus) && ~ismember(mainfluxnames{1}{i2},tracked)
                        track_clus=[track_clus mainfluxnames{1}{i2}];
                    end
                end
            end
        end
    end
    tracked=[tracked fluxes(i1).start];
end

while ~isempty(track_clus)
    track_clus_new={};
    for i1=1:length(track_clus)
        if check_cluster(track_clus{i1}) && ~ismember(track_clus{i1},tracked)
            start_ch=solve_charge(track_clus{i1});
            if start_ch==ch
                if exist('IPR','var')
                    [~,~,mainfluxnames,mainfluxes,~]=plotflux(name_mon,C_mon,track_clus{i1},'to',0,...
                        'C_mon_notation',C_mon_notation,'crit',crit_clust,'dp',dirpath,'suffix',suffix,...
                        'lsame_ch','IPR',IPR);
                else
                    [~,~,mainfluxnames,mainfluxes,~]=plotflux(name_mon,C_mon,track_clus{i1},'to',0,...
                        'C_mon_notation',C_mon_notation,'crit',crit_clust,'dp',dirpath,'suffix',suffix,...
                        'lsame_ch');
                end
                for i2=1:length(mainfluxes{1})
                    if mainfluxes{1}(i2)>=thflux && ~strcmpi(mainfluxnames{1}{i2},'others')
                        k=k+1;
                        fluxes(k).end=track_clus{i1};
                        fluxes(k).start=mainfluxnames{1}{i2};
                        fluxes(k).val=mainfluxes{1}(i2);
                        if ~ismember(mainfluxnames{1}{i2},track_clus_new) && ~ismember(mainfluxnames{1}{i2},tracked)
                            track_clus_new=[track_clus_new mainfluxnames{1}{i2}];
                        end
                    end
                end
            end
        end
        tracked=[tracked track_clus{i1}];
    end
    track_clus=track_clus_new;
end

figure(1)
clf(1)
set(gcf,'Position',[50 500 500 420]);
set(gca,'Position',[0.1 0.12 0.8 0.8]);
set(gca,'FontSize',12,'FontWeight','normal');
set(gcf,'Color','white')
box on
hold on
xlabel(['Molecules ',acidname])
ylabel(['Molecules ',basename])
title(titlestr,'FontWeight','normal')
if maxacids>maxbases
    maxlim=maxacids;
else
    maxlim=maxbases;
end
xlim([0 2*maxlim]);
ylim([0 2*maxlim]);
line([0 maxacids+1],[maxbases+1 maxbases+1],'Color','k','LineWidth',2.5,'LineStyle','--');
line([maxacids+1 maxacids+1],[0 maxbases+1],'Color','k','LineWidth',2.5,'LineStyle','--');
% Set the arrow colors and widths
arrowfluxes=[fluxes.val];
for i=1:length(arrowfluxes)
    if ~check_cluster(fluxes(i).start) || ~check_cluster(fluxes(i).end)
        arrowfluxes(i)=nan;
    end
end
step=(log10(max(arrowfluxes,[],'omitnan'))-log10(min(arrowfluxes,[],'omitnan')))/nlevels;
limits=log10(min(arrowfluxes,[],'omitnan')):step:log10(max(arrowfluxes,[],'omitnan'));
Lwidth=1:0.5:0.5*(nlevels+1);
colors=makeColorMap([0 0.45 0.74],[1 0 0],nlevels); % blue -> red
cmap=makeColorMap([0 0.45 0.74],[1 0 0]);

maxacidsend=0;
maxbasesend=0;
for i=1:size(fluxes,2)
    if check_cluster(fluxes(i).start)
        start_ch=solve_charge(fluxes(i).start);
        end_ch=solve_charge(fluxes(i).end);
        if start_ch~=ch && end_ch==ch
            % Just mark the source to the figure, if it's not of the same
            % charge (unless it's just a neutral monomer)
            if ~strcmpi(fluxes(i).start,['1',acidname]) && ~strcmpi(fluxes(i).start,['1',basename])
                [nacids,nbases]=countnumbers(fluxes(i).end,acidname,basename,acidname_neg);
                if start_ch==1
                    chcolor='k';
                elseif start_ch==2
                    chcolor='b';
                else
                    chcolor='r';
                end
                text(nacids+0.5,nbases+0.5,fluxes(i).start,'EdgeColor',chcolor,'Fontsize',16,'Linewidth',2);
            end
        elseif start_ch==ch && end_ch==ch
            % Draw an arrow
            [nacids,nbases]=countnumbers(fluxes(i).start,acidname,basename,acidname_neg);
            startcoord=[nacids+0.5,nbases+0.5];
            [nacids,nbases]=countnumbers(fluxes(i).end,acidname,basename,acidname_neg);
            endcoord=[nacids+0.5,nbases+0.5];
            if nacids>maxacidsend, maxacidsend=nacids; end
            if nbases>maxbasesend, maxbasesend=nbases; end
            nlev=nan;
            for j=1:length(limits)-1
                if log10(fluxes(i).val)>=limits(j) && log10(fluxes(i).val)<=limits(j+1)
                    nlev=j;
                    break
                end
            end
            arrow(startcoord,endcoord,'Linewidth',Lwidth(nlev),'Edgecolor',colors(nlev,:),'Facecolor',colors(nlev,:))
        end
    else
        % Mark the source to the figure, if it's not a cluster
        if ~strcmpi(fluxes(i).start,'source')
            [nacids,nbases]=countnumbers(fluxes(i).end,acidname,basename,acidname_neg);
            text(nacids+0.5,nbases+0.5,fluxes(i).start,'EdgeColor','m','Fontsize',16,'Linewidth',2);
        end
    end
end
pause
if maxacidsend>maxbasesend
    maxlim=maxacidsend;
else
    maxlim=maxbasesend;
end
xlim([0 maxlim+1]);
ylim([0 maxlim+1]);
colormap(cmap);
for i=1:length(limits)
    cbtick{i}=get_es_str(10^limits(i),0);
end
cb=colorbar('YTick',0:1/nlevels:1,'YTickLabel',cbtick,'FontWeight','normal','FontSize',12);
ylabel(cb,'Flux (cm^{-3}s^{-1})','FontWeight','normal','FontSize',12)
set(gca,'XTick',(0:maxlim)+0.5,'XTickLabel',0:maxlim);
set(gca,'YTick',(0:maxlim)+0.5,'YTickLabel',0:maxlim);

% Grid lines
g_x=get(gca,'XTick')+0.5;
g_x=[min(g_x-1),g_x];
g_y=get(gca,'YTick')+0.5;
g_y=[min(g_y-1),g_y];
for i=1:length(g_x)
   h=plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'Color',[0.38 0.38 0.38],'LineWidth',1,'LineStyle','-'); % y grid lines
   uistack(h,'bottom')
end
for i=1:length(g_y)
   h=plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'Color',[0.38 0.38 0.38],'LineWidth',1,'LineStyle','-'); % x grid lines
   uistack(h,'bottom')
end

hold off

end

function [nacids,nbases]=countnumbers(namelabel,acidname,basename,acidname_neg)
% Find the numbers of acids and bases for the matrix plot
nacids=0;
nbases=0;
[molnames,nmols]=parse_cluster(namelabel);
if ismember(acidname,molnames)
    nA=strcmp(molnames,acidname);
    nacids=nmols(nA);
end
if ismember(basename,molnames)
    nB=strcmp(molnames,basename);
    nbases=nmols(nB);
end
% Add the ion to the number of acids, if this is a negative cluster
if ismember(acidname_neg,molnames)
    nA_neg=strcmp(molnames,acidname_neg);
    if nmols(nA_neg)~=1
        error(['Number of negative ions is ',num2str(nmols(nA_neg)),'??'])
    end
    nacids=nacids+nmols(nA_neg);
end

end