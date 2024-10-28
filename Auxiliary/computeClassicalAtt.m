function alpha=computeClassicalAtt(model,c,f)
%alpha=computeClassicalAtt(model,c,f)
%A function that compute the coeff of classical attenuation alpha at the
%frequency f(Hz) for the gas and P,T conditions contained in model for a
%speed of sound c(m/s);
 
try
gas_p=model.gasT;
catch
    gas_p=readtable("gas_parametersSIModel.csv");
end
molmass=gas_p.molar_mass'*gas_p.MolarFract;
alpha=2*pi^2*f.^2./(model.rho*c.^3)*(4/3*model.etha+(model.cp/model.cv-1)*model.kappa/model.cp*molmass);



end