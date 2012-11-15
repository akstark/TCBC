function [y, t]=kinetics_solver(T,P,time)
%This function integrates a cantera kinetic mechanism at the prescribed
%Temperature, T, Pressure, P, and time, time as entered above (T,P,time)
%returning a matrix of N rows and M columns, N being the time resolution 
%used to define tspan below, and N being the total number of species in 
%the cantera gas object.
%The solver used is ode15s which calls the function contained below f, to
%return the net production rates of the species present in the gas mixture.
%
%Created by Addison Killean Stark on 9/16/2012



cleanup;                                                                   
% clean up cantera
                                   
                                 
tspan = linspace(0,time,100)';                                             
% divide the time of run for plot             

gas = IdealGasMix('pyrolysis.xml');                                        
% create ideal gas object from cantera cti

y0=zeros(1,nTotalSpecies(gas));                                            
% initialize vector of mass fractions to zeros
y0(324)=1;                                                                 
% add in species of interest... THIS COULD BE DONE BY NAME
y0=y0/sum(y0);
options=odeset('RelTol',1.e-3,'AbsTol',1.e-10,'Stats','on');               
% options for ode solver
[t,y]=ode15s(@f,tspan,y0,options,gas,T,P);                                 
% ode solver ode15s(@function ODE vector, tspan from above, initial condition, options from above, gas item from above)
plot(t,y)                                                                  

end

function dydt = f(t,y,gas,T,P)                                             
% the function f takes the arguments f(time,current values of y, gas) and returns dydt
                                                                           
% the temperature in Celsius.
Temp=273+T+t*0;                                                            
% set temperature in K, this needs to be made into dT/dt function
Press=P;                                                                   
% set pressure in atm, this needs to be made into dP/dt function
set(gas,'T',Temp,'P',Press*oneatm,'Y',y);                                  
% given T and P and y this creates an updated gas at t, all kinetic information is held here
dydt = netProdRates(gas).*molarMasses(gas);                                                  
% net production rates (kmol/(m^3 s)) are returned as the evolution of mass fractions 
end

