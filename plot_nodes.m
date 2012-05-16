function plot_nodes(XZ)
%PLOT_NODES(XZ)
%  XZ:  matrix of xz-coordinates

p = plot(XZ(:,1),XZ(:,2),'-o');

% axes formatting
axis equal;

% plot formatting
set(p, 'Color', 'black', ...
       'LineWidth', 2);

