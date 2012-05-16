function plot_mode_shape(nodes_xz, nodes_m_xz, i, f)
%PLOT_MODE_SHAPE(NODES_XZ, NODES_M_XZ)
%  NODES_XZ:   matrix of xz-coordinates of undisplaced nodes
%  NODES_M_XZ: matrix of xz-coordinates of discplaced nodes
%  I:          mode nr.
%  F:          mode frequency

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

