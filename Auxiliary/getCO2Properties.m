function [cv,cp,eta,kappa]=getCO2Properties(T)
%A function that return some fluid properties of co2 at the temperature T in K based on
%extrapolation from webbook nist data 

cv=0.0486*T+14.3683;
cp=0.0485*T+22.6926;
kappa=7.29e-5*T-0.0052;
eta=4.8468e-8*T+4.605e-7;

end

