function [ nodes_m ] = get_mode_shape( nodes, V )
%GET_MODE_SHAPE Returns displaced node coordinates of a vibration mode
%using the corresponding eigenvector
%   Input arguments:  nodes   - xz-coordinates of all nodes
%                     V       - eigenvector of modal displacements
%   Output arguments: nodes_m - xz-coordinates of displaced nodes

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

