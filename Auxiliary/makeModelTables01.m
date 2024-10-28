function makeModelTables01(model)
%A function that create and save a tables containing all the molecule
%selected in the model, their vibration modes and their molar fraction.

%Loading parameters
load_physic_constant
gas_p=readtable("gas_parametersSI.csv");
vib_p=readtable("vibration_mode_parametersSI.csv");

gas_p_mod=gas_p(ismember(gas_p.Name,model.listOfMole),:);
vib_p_mod=vib_p(ismember(vib_p.MolName,model.listOfMole)&ismember(vib_p.VibrationMode,model.listOfMode),:);


[ia ib]=ismember(gas_p_mod.Name,unique(model.listOfMole,'stable'));
gas_p_mod.MolarFract=model.MolarFrac(ib)';
gas_p_mod.N(:)=model.P/(R*model.T)*Na;%*model.MolarFrac;
[ia ib]=ismember(vib_p_mod.MolName,gas_p_mod.Name);
vib_p_mod.MolarFract=gas_p_mod.MolarFract(ib);

%Make table for gas and vibration mode corresponding to the order of the
%model
model_tab_gas=table(cellstr(unique(model.listOfMole','stable')));
model_tab_gas.Properties.VariableNames=["Name"];

model_tab_vib=table(cellstr(model.listOfMole'),cellstr(model.listOfMode'));
model_tab_vib.Properties.VariableNames=["MolName","VibrationMode"];

gas_p_mod=join(model_tab_gas,gas_p_mod,'Key',"Name");

n=length(model.listOfMode);
vib=vib_p_mod;
for i=1:n
    line=vib_p_mod(string(vib_p_mod.MolName)==model.listOfMole(i)&string(vib_p_mod.VibrationMode)==model.listOfMode(i),:);
    vib(i,:)=line;
end
writetable(vib,"Tables\vibration_mode_parametersSIModel.csv");
writetable(gas_p_mod,"Tables\gas_parametersSIModel.csv");

end