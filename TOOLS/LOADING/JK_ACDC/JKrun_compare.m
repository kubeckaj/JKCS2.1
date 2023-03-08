%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, input2.inp dG.txt, dG2.txt
%2) in ACDC folder, adjust: ACDCinit, ACDCsource, ACDCinit2, ACDCsource2
%3) run 
%4) plot some properties: plotJT, plotCT for the last ACDCrun
%5) use stop to interupt after first run if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
Temp=298.15; %Temperature

rh=100;

ACDCinit(Temp);                      %adjust ACDCinit
ACDCsource(10^12,rh2cmmc(rh,Temp))   %adjust ACDCsource
[T,C,J_out]=ACDCrun(); 
J_out(end)

%stop  %uncomment  stop in between if you need
 
ACDCinit2(Temp,rh);                  %adjust ACDCinit
ACDCsource2(10^12);                  %adjust ACDCsource
[T,C,J_out]=ACDCrun(); 
J_out(end)

%get_evap; get_coll; get_wl; get_cs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'