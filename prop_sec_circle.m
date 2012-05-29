function [ rhoA,EA,EI ] = calc_prop_circle(rho,E,mue,length,diameter)
%Calculation of Material Properties of a circle cross-section beam
%with constant cross-section
%   input Parameters: 
    % rho = material density (kg/m^3)
    % E = Youngs Modulus (MPa)
    % mue =  poissons constant (-)
    % diameter, height = beam dimesnions (m)

% cross-section area (m^2)
A=pi*(diameter^2)/4;    
% volume (m³)
V=A*length;
% mass (kg)
m=rho*V;
% line mass (kg/m)
rhoA = rho*A; 
% cross-section inertia (m^4)
I=pi*(diameter^4)/64;                     
% axial rigidity (N)
EA = E*10^6*A;                      
% bending rigidity (N*m^2)
EI = (E*10^6)/I;    

end

