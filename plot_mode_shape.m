function plot_mode_shape( nodes_xz, nodes_m_xz, i, f )
%PLOT_MODE_SHAPE - Plots the deformed and undeformed eigenmode-system
% This function plots the system in the undeformed and deformed state
% for a specified eigenmode. The deformed node coordinates plus the mode
% number and the mode eigenfrequency have to be specified.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_mode_shape.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : plot_mode_shape( nodes_xz, nodes_m_xz, i, f )
%
% Inputs :
%    nodes_xz   - x- and z-coordinates of the nodes of the undeformed
%                 system
%    nodes_m_xz - x- and z-coordinates of the nodes of the deformed
%                 system
%    i          - mode number
%    f          - mode eigenfrequency
%
% Outputs :
%    none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-17 12:57 CEST
% Last Modified : 2012-05-17 13:01 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%p = plot(nodes_xz(:,1),nodes_xz(:,2),'--x', ...
%         nodes_m_xz(:,1),nodes_m_xz(:,2),'-o');

% axes formatting
axis equal;
axis([0 2 -0.25 0.25]);

% title
title(sprintf('Mode %i: f = %6.1f Hz', [i f]));

% plot formatting
%set('Color', 'black', ...
      % 'LineWidth', 2);

line(nodes_xz(:,1), nodes_xz(:,2), ...
     'Color', [.5 .5 .5], ...
     'LineStyle', '--', ...
     'Marker', 's');

line(nodes_m_xz(:,1), nodes_m_xz(:,2), ...
     'Color', 'black', ...
     'LineWidth', 2, ...
     'Marker', 'o');

