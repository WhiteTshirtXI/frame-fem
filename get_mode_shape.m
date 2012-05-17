function [ nodes_m ] = get_mode_shape( nodes, V )
%GET_MODE_SHAPE - Return frame deformation for a specified eigenmode
% This function returns the x- and z-coordinates of all system nodes
% when deformed by a eigenmode specified through a eigenvector.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_mode_shape.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : [ nodes_m ] = get_mode_shape( nodes, V )
%
% Inputs :
%    nodes - x- and z-coordinates of all system nodes
%    V     - eigenvector of nodal displacements
%
% Outputs :
%    nodes_m - x- and z-coordinates of the displaced system nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-17 12:53 CEST
% Last Modified : 2012-05-17 12:56 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the eigenvector is a vector of size:
% (number of nodes)*(3 degrees of freedom)
% -> to get the displacements in x and z direction, every third row has
% to be skipped!

% x displacement
m_displacement_x = V(1:3:end);
% z displacement
m_displacement_z = V(2:3:end);

% assemble
m_displacement_xz = [m_displacement_x m_displacement_z];

% normalize with maximum value
m_displacement_xz = 0.1*m_displacement_xz ./ max(m_displacement_xz(:));

% superposition
nodes_m = nodes + m_displacement_xz;  


end

