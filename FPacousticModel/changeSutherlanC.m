function oldC=changeSutherlanC(molname,C)
%A function that replace the Sutherland constant for the species
%"molname" by the value C.
%% Check the current directory
list=ls;
list=deblank(string(list));
flag=sum(list=="PhysicalPropertiesTables");
if flag==0
    %disp("Error: the function must be executed in a directory that contains a folder named Tables with all the gas and vibration mode tables in it")
    %return
    cd("D:\m.gillier\Documents\MATLAB\AcousticPropagationModel\Attenuation\Data")
end

load_physic_constant
gas_p=readtable("gas_parameters.csv");
line=gas_p(gas_p.Name==molname,:);
if height(line)==1
    oldC=gas_p(gas_p.Name==molname,:).SutherlandC;
        gas_p(gas_p.Name==molname,:).SutherlandC=C;
    
    writetable(gas_p,"PhysicalPropertiesTables\gas_parameters.csv");

else
    disp("Error: the selected molecule is not present in the current model.")
    return
end
%% Converting to SI
gas_parameter=readtable("gas_parameters.csv");
vibration_mode_parameter=readtable("vibration_mode_parameters.csv");
gp_si=convertTable(gas_parameter);
vmp_si=convertTable(vibration_mode_parameter);
writetable(vmp_si,"PhysicalPropertiesTables\vibration_mode_parametersSI.csv");
writetable(gp_si,"PhysicalPropertiesTables\gas_parametersSI.csv");

end