function newt=convertTable(table)
%a function that convert a table of gas molecule parameters or vibration mode parameter
%in IS units.
load_physic_constant
newt=table;
if string(table.Properties.VariableNames{1})=="MolName"%the table is the vibration mode one
    newt.nu=newt.nu*100;%cm-1 to m-1
    newt.VibAmplCoeff=newt.VibAmplCoeff*1000*Na;%from amu-1((g/mol)-1) to kg???????;
elseif table.Properties.VariableNames{1}=='Name'%the table is the gas parameters one
    newt.molar_mass=newt.molar_mass/1000;%from g/mol to kg/mol
    newt.sigma=newt.sigma*1e-10;%from Angstrom to meters
    newt.epsilon=newt.epsilon*kb;%in J ,from  epsilon/kb to epsilon
end


end