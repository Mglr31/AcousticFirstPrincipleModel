function changeLJparameter(molname,sigma,epsilon)
%A function that replace the Lennard-Jones parameters for the species
%"molname" by the value of sigma(Angstrom) and epsilon(e/KbT, K) in the
%table gas_parametersSI.If sigma(or epsilon) equal 0 the value will
%not be changed.
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
sigma=sigma;
epsilon=epsilon;
gas_p=readtable("gas_parameters.csv");
line=gas_p(gas_p.Name==molname,:);
if height(line)==1
    if sigma~=0
        gas_p(gas_p.Name==molname,:).sigma=sigma;
    end
    if epsilon~=0
        gas_p(gas_p.Name==molname,:).epsilon=epsilon;
    end
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