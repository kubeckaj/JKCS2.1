function y = pa2cmmc(p,T) 
  % Avogadro's constant
  Na=6.02214179e23;
  % Boltzmann constant [J/K]
  kB=1.3806504e-23;
  R  = 8.314;
  
  y = 10^-6*Na*p/R./T;
end