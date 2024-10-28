function newmode=addDerivativePArameter(model)
%A function that add to a  model some parameetrs computed from the model
%input such as the density ,cp,cv,cvib...
load_physic_constant
%mol_p=readtable("gas_parametersSIModel.csv");
mol_p=model.gasT;

MolMasses=mol_p.molar_mass;
model.rho=sum(model.P/(model.T*R)*MolMasses.*mol_p.MolarFract);
model.MolarMass=MolMasses'*mol_p.MolarFract;
[~, lib]=ismember(unique(model.listOfMole,'stable'),mol_p.Name);%getting the list of thermal cond and viscosity for the model
ethas_model=mol_p.etha(lib);
kappas_model=mol_p.kappa(lib);
molMasse_model=mol_p.molar_mass(lib);
%Computing the mixture thermal cond and viscosity
n=height(mol_p);
phi=zeros(n);
for a=1:n
    for b=1:n
        phi(a,b)=(1+sqrt(ethas_model(a)/ethas_model(b))*(molMasse_model(b)/molMasse_model(a))^0.25)^2;
        phi(a,b)=phi(a,b)/sqrt(8*(1+molMasse_model(a)/molMasse_model(b)));
    end
end

model.etha=wilkeExpression(model.MolarFrac,ethas_model,phi);
model.kappa=wilkeExpression(model.MolarFrac,kappas_model,phi);




%vib_p=readtable("vibration_mode_parametersSIModel.csv");
vib_p=model.vibT;
[~, lib]=ismember(model.listOfMole,vib_p.MolName);
model.MolarFractMode=vib_p.MolarFract(lib)';
model.cvib=arrayfun(@(x,y) specificHeatVib(x,y,model.T,model),model.listOfMole,model.listOfMode);
model.cv=sum(mol_p.cv.*mol_p.MolarFract);
model.cp=sum(mol_p.cp.*mol_p.MolarFract);
model.c=sqrt(1.4*R*model.T);


newmode=model;
end