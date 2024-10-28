function [cv,cp,eta,kappa]=getH2OProperties(T)
%A function that return some fluid properties of H20 at the temperature T in K based on
%extrapolation from webbook nist data 

cv=25.519;
cp=33.908;
kappa=0.016815;
eta=8.9700e-06;
end