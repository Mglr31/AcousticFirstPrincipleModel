function a=getThermalDiffusivity(model)
% a=getThermalDiffusivity(model)

cp=model.cp/model.MolarMass;%cp in J/kg*K
a=model.kappa/model.rho/cp;


end