Lwidth=100; %Box Size
Fsize=12;   %Font Size
figure();   %Initiate new figure

%%% SELECT DATA TO  PLOT
%monY_coll_mat monX_coll_mat evap_tot_mat deltag_mat_new deltag_mat
data=flip(monX_coll_mat);

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
colormap('parula');brighten(0.7)

%%% ADD LEGEND
cb=colorbar('FontSize',Fsize,'FontWeight','demi');
data_label=['{\it c}_{' Xname '}{\it\beta}_{' Xname '} [s^{-1}]'];
ylabel(cb,data_label,'LineWidth',Lwidth,'FontSize',Fsize,'FontWeight','demi')
ylabel(["\fontsize{20}" Yname])
xlabel(["\fontsize{20}" Xname])
title("\fontsize{20}Collision rate with "+Xname)

%%% SETUP FIGURE SIZE AND POSITION
figpos=[0 0 (max_X+1)*Lwidth (max_Y+1)*Lwidth];set(gcf,'Position',figpos);
set(gca,'XTick',(1:max_X+1)+0.5,'XTickLabel',0:max_X,'FontSize',Fsize,'FontWeight','demi');
set(gca,'YTick',(1:max_Y+1)+0.5,'YTickLabel',0:max_Y,'FontSize',Fsize,'FontWeight','demi');
set(gca,'Position',[0.08 0.18 0.95*0.78 0.95*0.75]);
pos=get(gca,'Position');

%%% INSERT BOX LABELS
width=pos(3)/(max_X+1);
height =pos(4)/(max_Y+1);
for i=1:max_X+1
  for j=1:max_Y+1
      if isnan(data(j,i))
          str="";
      else
          %if max(max(log10(data)))<2
          %  str=sprintf('%0.2f',data(j,i));
          %else
            str=get_es_str(10^log10(data(j,i)),0);
          %end
      end
      annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
                            'string',str,'LineStyle','none','HorizontalAlignment', ...
                            'center','VerticalAlignment','middle','FontWeight','demi','FontSize',Fsize);
  end
end