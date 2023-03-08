cd RUN
source=zeros([size(C(end,:),2),1]);
[flux_2d, coll_evap_2d, sources, flux_3d, flux_out, clust_flux]=dofluxes(C(end,:)*10^6,K,E,WL,CS,source,0,0,0);

Lwidth=100; %Box Size
Fsize=12;   %Font Size
figure();   %Initiate new figure
dif=10;     %Amplitude difference: 10^(MAX-dif) is visible
mult=6;     %Arrow max size

%%% SELECT DATA TO  PLOT
%monY_coll_mat monX_coll_mat evap_tot_mat deltag_mat_new deltag_mat
data=flip((monY_coll_mat+monX_coll_mat)./evap_tot_mat);

%%% ADD 1 row and 1 column
dataP=nan(size(data)+[1,1]);
dataP(1:size(data,1),1:size(data,2))=data;
data=dataP;

%%% FORMAT OF DATA
%if max(max(log10(data)))<2
%  pcolor_mat=[data NaN(size(data,1),1); NaN(1,size(data,2)+1)];
%else
  pcolor_mat=[log10(data) NaN(size(data,1),1); NaN(1,size(data,2)+1)];
%end

%%% COLOR MAP
pcolor(pcolor_mat);

%%% TYPE OF COLOR
%colormap('cool');
colormap('parula');brighten(0.7);

%%% ADD LEGEND
cb=colorbar('FontSize',Fsize,'FontWeight','demi');
data_label=['({\it c}_{' Xname '}{\it\beta}_{' Xname '}+{\it c}_{' Yname '}{\it\beta}_{' Yname '}) / \Sigma{\it\gamma} [  ]'];
ylabel(cb,data_label,'LineWidth',Lwidth,'FontSize',Fsize,'FontWeight','demi')
ylabel(["\fontsize{20}" Yname])
xlabel(["\fontsize{20}" Xname])
title("SYSTEM FLUXES")

%%% SETUP FIGURE SIZE AND POSITION
figpos=[0 0 (max_X+2)*Lwidth (max_Y+2)*Lwidth];set(gcf,'Position',figpos);
set(gca,'XTick',(1:max_X+2)+0.5,'XTickLabel',0:max_X+1,'FontSize',Fsize,'FontWeight','demi');
set(gca,'YTick',(1:max_Y+2)+0.5,'YTickLabel',0:max_Y+1,'FontSize',Fsize,'FontWeight','demi');
set(gca,'Position',[0.08 0.18 0.95*0.78 0.95*0.75]);
pos=get(gca,'Position');

%%% INSERT BOX LABELS
width=pos(3)/(max_X+2);
height =pos(4)/(max_Y+2);
maxflux=max(max(log10(flux_2d)));
mapf=flip(map);
for i=1:max_X+1
  for j=1:max_Y+1
      %arrow to up
      kS1=NaN;
      kS2=NaN;
      for k=1:length(clust_acdc)
         if compare_clusters(mapf{j,i},clust_acdc{k})
             kS1=k;
         end
         if j+1<max_Y+2
           if compare_clusters(mapf{j+1,i},clust_acdc{k})
             kS2=k;
           end
         end
      end
      if isnan(kS2)
          kS2=size(flux_2d,1)-1;
      end
      if isnan(kS1)
          fW=-Inf;
      else
          diff=flux_2d(kS1,kS2)-flux_2d(kS2,kS1);
          fW=(log10(abs(diff))-maxflux+dif)/dif;
      end
      if fW<0
         fW=NaN;
      end
      if not(isnan(fW)) && not(isinf(fW))
         ar = annotation('arrow');
         if diff>0
           ar.X=[pos(1)+width*(i-1+0.5),pos(1)+width*(i-1+0.5)];
           ar.Y=[pos(2)+height*(j-1+0.5),pos(2)+height*(j+0.5)];
         else
           ar.X=[pos(1)+width*(i-1+0.5),pos(1)+width*(i-1+0.5)];
           ar.Y=[pos(2)+height*(j+0.5),pos(2)+height*(j-1+0.5)];
         end
         ar.LineWidth=mult*fW;
      end
      %arrow to right
      kS1=NaN;
      kS2=NaN;
      for k=1:length(clust_acdc)
         if compare_clusters(mapf{j,i},clust_acdc{k})
             kS1=k;
         end
         if i+1<max_X+2
           if compare_clusters(mapf{j,i+1},clust_acdc{k})
             kS2=k;
           end
         end
      end
      if isnan(kS2)
          kS2=size(flux_2d,1)-1;
      end
      if isnan(kS1)
          fW=-Inf;
      else
          diff=flux_2d(kS1,kS2)-flux_2d(kS2,kS1);
          fW=(log10(abs(diff))-maxflux+dif)/dif;
      end
      if fW<0
         fW=NaN;
      end
      if not(isnan(fW)) && not(isinf(fW))
         ar = annotation('arrow');
         if diff>0
           ar.X=[pos(1)+width*(i-1+0.5),pos(1)+width*(i+0.5)];
           ar.Y=[pos(2)+height*(j-1+0.5),pos(2)+height*(j-1+0.5)];
         else
           ar.X=[pos(1)+width*(i+0.5),pos(1)+width*(i-1+0.5)];
           ar.Y=[pos(2)+height*(j-1+0.5),pos(2)+height*(j-1+0.5)];
         end
         ar.LineWidth=mult*fW;
      end
  end
end
cd ..