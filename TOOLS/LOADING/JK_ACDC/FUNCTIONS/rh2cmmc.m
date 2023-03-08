function y = rh2cmmc(rh,temp) 
  % Avogadro's constant
  Na=6.02214179e23;
  % Boltzmann constant [J/K]
  kB=1.3806504e-23;
  R  = 8.314;
  
  % Water: RH
  y=rh/100*exp(-2991.2729*temp^(-2)-6017.0128/temp+18.87643854-0.028354721*temp...
    +0.17838301e-4*temp^2-0.84150417e-9*temp^3+0.44412543e-12*temp^4+2.858487*log(temp));
  y=10^-6*Na*y/R./temp;
end