% A MATLAB script that shows how to use the attenuationModel.m function to compute the attenuation and speed of sound 
% in a gase mixture, in this example the standard Martian atmosphere 

%% Making an atmosphere
%Example : standard dry Martian atmosphere
atmMars.listOfMole=["N2","CO2","CO2","CO2","Ar","O2"];%List of molecules making the Martian atm, repeated when multiple 
% vibration modes are considered
atmMars.listOfMode=["nu", "nu1", "nu2","nu3","none","nu"];%Name of the vibration modes corresponding to the above molecules
atmMars.MolarFrac=[0.0270 0.9500 0.0160 0.0013];%Molar fraction of each species of the atmosphere ( "N2","CO2","Ar","O2")
atmMars.T=240;%Temperature in K
atmMars.P=740;%Pressure in Pa

%% Computing alpha and c for all frequencies
%It is best to have an Internet connection so that the thermophysical
%properties of each gas at the required P and T can be pulled from the NIST online database.
model=attenuationModel(atmMars);

%% Plotting the result
f1=figure();
%Plotting the attenuation coefficient
loglog(model.f,model.alpha,'DisplayName',"\alpha",'LineWidth',2)%Total attenuation
hold on
loglog(model.f,model.alpha_c,'DisplayName',"\alpha_{class}")%Classical attenuation
loglog(model.f,model.alpha_r,'DisplayName',"\alpha_{mol}")%Molecular (or vibrational) attenuation

legend('Location','northwest')
grid on
xlabel("Frequency [Hz]")
ylabel("Attenuation coefficient [m^{-1}]")
xlim([0.1 1e6])
set(gca,'FontSize',20)

f2=figure();
%Plotting the speed of sound 
loglog(model.f,model.c,'LineWidth',2)

grid on
xlabel("Frequency [Hz]")
ylabel("Speed of sound [m.s^{-1}]")
xlim([0.1 1e6])
set(gca,'FontSize',20)
