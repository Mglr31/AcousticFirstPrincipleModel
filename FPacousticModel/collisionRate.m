function z=collisionRate(nmole1,nmole2,T)
%z=collisionRate(nmole1,nmole2,T)
%Compute the collision rate z (s-1) between the species nmole1 and
%nmole2(nmole="H20","N2","CH4",...)
%at the temperature T(K).

%% Loading parameters
load_physic_constant

if isstring(nmole1)
    disp(strcat("collisionRate: Reading table for",nmole1))

    gas_p=readtable("gas_parametersSIModel.csv");

    mol1=gas_p(gas_p.Name==nmole1,:);

else
    mol1=nmole1;
end
if isstring(nmole2)
    disp(strcat("collisionRate: Reading table for",nmole2))

    gas_p=readtable("gas_parametersSIModel.csv");


    mol2=gas_p(gas_p.Name==nmole2,:);
else
    mol2=nmole2;
end
%% Computing Z

z=2*mol2.N*((mol1.sigma+mol2.sigma)/2)^2*sqrt(2*pi*kb*T*(mol1.mass+mol2.mass)/(mol1.mass*mol2.mass));

end