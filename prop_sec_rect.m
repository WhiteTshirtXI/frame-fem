function [ rhoA,EA,EI ] = prop_sec_rect(rho,E,mue,length,width,heigth)
%Calculation of Material Properties of a rectangular cross-section beam
%with constant cross-section
%   input Parameters: 
    % rho = material density (kg/m^3)
    % E = Youngs Modulus (MPa)
    % mue =  poissons constant (-)
    % length, width, height = beam dimesnions (m)

% volume (m³)
V=heigth*width*length;
% mass (kg)
m=rho*V;
% cross-section area (m^2)
A=heigth*width;
% line mass (kg/m)
rhoA = rho*A; 
% cross-section inertia (m^4)
I=(heigth^3)*width/12;                     
% axial rigidity (N)
EA = E*10^6*A;                      
% bending rigidity of a beam (N*m^2)
EI = E*10^6 * I;
% bending rigidity of a plate (N*m^2)
%EI = (E*10^6)/(1-mue^2)*I;

end

