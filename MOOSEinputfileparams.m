%this file calculated input parameters for MOOSE 

% DCPD PDMS Wire  -- order of entries



rho = [1 2 3];

k = [2 3 4];

Cp = [3 4 6];

Hr = 2;
%% 

L = 3; %Characterictic length

A = 1.91e4;  %from kessler's paper

Ttrig = 339.44;  %Ignition temperatue in degree C

Tinit = 20;  %Intial Temperature in degree C

%% 

%DCPD-PDMS-WIRE
T_condutivity = k./rho(1)/Cp(1)/A/L^2;

Hrinput = Hr/Cp(1)/(Ttrig-Tinit);

rhocp = rho.*Cp./(rho(1)*Cp(1));

%% Print



