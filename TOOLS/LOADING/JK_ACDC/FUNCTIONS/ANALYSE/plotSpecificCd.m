%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ ACDC CLUSTERS
fp=fopen('../RUN/driver_acdc.m','r');
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

newC=[];
newN=[];
for j=1:size(C,2)
    if C(end,j) > 1000
       newC=[newC,C(end,j)]; 
       newN=[newN,clust_acdc(j)];
    end
end
newN=["(SA)_1","(MSA)_1","(DMA)_1","(SA)_1(DMA)_1","(SA)_2(DMA)_1","(SA)_2(DMA)_2","(SA)_3(DMA)_2","(SA)_3(DMA)_3","(SA)_4(DMA)_4","(MSA)_2(DMA)_2","(SA)_1(MSA)_1(DMA)_1","(SA)_1(MSA)_1(DMA)_2","(SA)_2(MSA)_1(DMA)_3","(SA)_3(MSA)_1(DMA)_3","  (SA)_3(MSA)_1(DMA)_4"]

figure()
bar(newC)
hold on
box on
set(gca,'YScale','log')
set(gcf,'Color','white')
xlim([0,length(newC)+1])
set(gca,'XTick',1:length(newN),'XTickLabel',newN);            
ax=gca;    
ax.XTickLabelRotation=90;
ylabel('{\itC}_{steady-state} [cm^{-3}]')