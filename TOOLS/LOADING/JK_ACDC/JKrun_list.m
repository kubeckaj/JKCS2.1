%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, dG.txt (or dHdS.txt)
%2) in ACDC folder, adjust: ACDCinit, ACDCsource
%3) run 
%4) plot VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
temp=280; %Temperature

results=[];
VARrange=6:0.2:8;
VARrange=10.^VARrange;
%VARrange=[0.01 0.1 1 10 100 1000 10000];
for VAR=VARrange
    fprintf("->"+num2str(VAR));
    
    sa_conc=VAR; %cmmc=cm^-3 
    am_conc=1000; %ppt
    
    ACDCinit(temp);                  %adjust ACDCinit
    ACDCsource(sa_conc,ppt2cmmc(am_conc,temp)) %adjust ACDCsource
    [T,C,J_out]=ACDCrun();  
    results=[results J_out(end)];
end

%PLOT
figure()
loglog(VARrange,results)
title("SA-W nucleation")
xlabel("SA concentration [cm^{-3}]")
ylabel("nucleation rate [cm^{-3}s^{-1}]")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'
            
