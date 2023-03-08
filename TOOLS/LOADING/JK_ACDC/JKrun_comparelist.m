%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, input2.inp dG.txt, dG2.txt
%2) in ACDC folder, adjust: ACDCinit, ACDCsource, ACDCinit2, ACDCsource2
%3) run 
%4) plot VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
Temp=298.15; %Temperature



results1=[];
results2=[];
VARrange=-2:0.5:2.5;
VARrange=10.^VARrange;
for VAR=VARrange
    fprintf("->"+num2str(VAR));
    rh=VAR;
    
    pause(2)
    ACDCinit(Temp);                      %adjust ACDCinit
    ACDCsource(10^12,rh2cmmc(rh,Temp))   %adjust ACDCsource
    [T,C,J_out]=ACDCrun(); 
    J_out(end);
    results1=[results1 J_out(end)];
    
    pause(2)
    ACDCinit2(Temp,rh);                  %adjust ACDCinit
    ACDCsource2(10^12);                  %adjust ACDCsource
    [T,C,J_out]=ACDCrun(); 
    J_out(end);
    results2=[results2 J_out(end)];
end

%PLOT
figure()
loglog(VARrange,results1,'ro-')
hold on
loglog(VARrange,results2,'b+-')
%xlim([1 1000]);
%ylim([10^14 10^21]);
title("Effect of humidity on SA nucleation")
xlabel("humidity [%]")
ylabel("nucleation rate [cm^{-3}s^{-1}]")

%get_evap; get_coll; get_wl; get_cs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'