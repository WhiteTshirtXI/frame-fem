function plot_nodes( XZ )
%PLOT_NODES - Plot the nodes and element of the frame structure
% This function uses the specified x- and z-coordinates of the system
% nodes to plot the frame structure.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_nodes.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : plot_nodes( XZ )
%
% Inputs :
%    XZ - Matrix of node xz-coordinates
%
% Outputs :
%    none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-17 12:43 CEST
% Last Modified : 2012-05-17 12:46 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = plot(XZ(:,1),XZ(:,2),'-o');

% axes formatting
axis equal;

% plot formatting
set(p, 'Color', 'black', ...
       'LineWidth', 2);

