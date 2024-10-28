
load_physic_constant
vibration_mode_parameter=readtable("vibration_mode_parametersSI.csv");

vibration_mode_parameter.omega=vibration_mode_parameter.nu*100*cl*2*pi;
vibration_mode_parameter.n=vibration_mode_parameter.nu*100*cl;%n is frequency in Hz(???)
vibration_mode_parameter.theta=h*vibration_mode_parameter.nu*100*cl/kb;%characteristic temperature for vibration(k)
vibration_mode_parameter.ID=strcat(vibration_mode_parameter.MolName,vibration_mode_parameter.VibrationMode);

writetable(vibration_mode_parameter,"Tables\vibration_mode_parameters.csv");



gas_parameter=readtable("gas_parameters.csv");
gas_parameter.mass=gas_parameter.molar_mass/Na/1000;%mass of a molecule in kg
writetable(gas_parameter,"Tables\gas_parameters.csv");
%% Converting to SI
gas_parameter=readtable("gas_parameters.csv");
vibration_mode_parameter=readtable("vibration_mode_parameters.csv");
gp_si=convertTable(gas_parameter);
vmp_si=convertTable(vibration_mode_parameter);
writetable(vmp_si,"vibration_mode_parametersSI.csv");
writetable(gp_si,"gas_parametersSI.csv");