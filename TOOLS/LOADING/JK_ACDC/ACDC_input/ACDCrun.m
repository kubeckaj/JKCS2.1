function [T,C,J_out] = ACDCrun()
%ACDCRUN Summary of this function goes here
%   Detailed explanation goes here

  time=10^-6;jump=10;convergence_test=0;
  %[C, T,converged, clust, Cf, J_out]=driver_acdc(time,'Sources_in','sources.txt');
  [C, T,converged, clust, Cf, ~, ~, J_out]=driver_acdc(time,'Sources_in','sources.txt');
  Jinit=J_out(end);
  while convergence_test==0
    time=time*jump;
    %[C, T,converged, clust, Cf, J_out]=driver_acdc(time,C,T,'Sources_in','sources.txt');
    [C, T,converged, clust, Cf, ~, ~, J_out]=driver_acdc(time,C,T,'Sources_in','sources.txt');
    Jend=J_out(end);
    if abs((Jend-Jinit)/Jinit) < 10^(-5) 
      convergence_test=1;
    end
    Jinit=Jend;
    %fprintf(""+num2str(Jend)+"_"+num2str(CHOMnow)+"_"+num2str(convergence_test)+"\n")
  end

end

