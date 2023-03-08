%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ ACDC CLUSTERS
fp=fopen('RUN/driver_acdc.m','r');
line=0;
while line ~= -1
    line=[fgetl(fp),' '];
    if ~isempty(strfind(line, 'clust = {'))
        line=strrep(line,'clust','clust_acdc');
        eval(line);
        break
    end
end
fclose(fp);
% -> clust_acdc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
bar(C(end,:))
hold on
box on
set(gca,'YScale','log')
set(gcf,'Color','white')
xlim([0,length(clust_acdc)+1])
%ylim([10^1,10^10])
set(gca,'XTick',1:length(clust_acdc),'XTickLabel',clust_acdc);            
ax=gca;    
ax.XTickLabelRotation=70;
ylabel('{\itC}_{steady-state} (cm^{-3})')
            