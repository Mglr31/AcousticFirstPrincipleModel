function f=plotLJfit(sigma,epsilon,rc,alpha)
r=linspace(0.01*sigma,1.5*sigma,3000);
Vlj=4*epsilon*((sigma*r.^-1).^12-(sigma*r.^-1).^6);
Estar=4*epsilon*((sigma*rc^-1)^12-(sigma*rc^-1)^6);

Vexp=(epsilon+Estar)*exp(-alpha*(r-rc))-epsilon;


f=figure

plot(r/sigma,Vlj/epsilon,'DisplayName',"Lennard-Jones Potential")
hold on
plot(r/sigma,Vexp/epsilon,'DisplayName',"Exponential Fit")
legend show
xline([1,rc/sigma])
xlabel("r/sigma")
ylabel("V/epsilon")
grid on

end