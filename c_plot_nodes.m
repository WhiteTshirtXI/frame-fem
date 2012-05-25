classdef c_plot_nodes < handle
%C_PLOT_NODES - Plot nodes of a fem-system
% This class aids in plotting all nodes and elements of a fem-system.
% The nodes can be drawn in an undeformed and (for example by eigenmodes
% or excitation forces) deformed kind of state, both simultaneously and
% seperately.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_plot_nodes.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Properties :
%    DEF_FIG_SET  - default figure settings
%    DEF_AXES_SET - default axes settings
%
% Methods :
%    c_plot_nodes  - constructor
%    plotDisplaced - plot displaced nodes for different displacement
%                    states
%    calcRowsCols  - calculate optimum row and column number of subplots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-23 14:51 CEST
% Last Modified : 2012-05-25 16:26 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)
        % default figure settings
        DEF_FIG_SET = struct('NumberTitle', 'Off');

        % default axes settings
        DEF_AXES_SET = struct('DataAspectRatio', [1 1 1]);

        % line styles
        LS_BLACK_SOLID = struct('Color', 'black', ...
                                'LineWidth', 2, ...
                                'Marker', 'o');
        LS_GREY_DASHED = struct('Color', [.4 .4 .4], ...
                                'LineStyle', '--', ...
                                'LineWidth', 1, ...
                                'Marker', 's');

        % default axis limits
        ax_limits = [-inf inf -inf inf];
    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_nodes - coordinates of system nodes
        function self = c_plot_nodes(p_nodes)

            % assign maximum and minimum coordinates to axis limits

            % minimum
            self.ax_limits([1 3]) = min(p_nodes) - [0.1 0.3];
            % maximum
            self.ax_limits([2 4]) = max(p_nodes) + [0.1 0.3];

        end

        % PLOT DISPLACED NODES FOR DIFFERENT DISPLACEMENT STATES
        %
        % Inputs:
        %   p_nodes  - undisplaced node coordinates [x1 z1 ;
        %                                            x2 z2 ;
        %                                             ...  ]
        %   p_dNodes - node displacement states [dx11 dz11 dx12 dz12 ;
        %                                        dx21 dz21 dx22 dz22 ;
        %                                                ...         ]
        function self = plotDisplaced(self, p_nodes, p_dNodes, p_Adj)

            % get number of displacement states
            nStates = size(p_dNodes,2)/2;

            % calculate optimum row and column numbers for subplot
            % display
            [n_rows, n_cols] = self.calcRowsCols(nStates);

            % create parent figure
            figure(self.DEF_FIG_SET);

            % for each displacement state, create a subplot
            for state = 1:nStates

                % get node displacements for current state
                displacements = p_dNodes(:,(2*state-1):(2*state));

                % calculate displaced nodes coordinates
                nodes_displaced = p_nodes + displacements;

                % create subplot
                subplot(n_rows, n_cols, state);
                hold on;
                axis equal;
                axis(self.ax_limits);

                % plot undisplaced and displaced nodes
                gplot(p_Adj, p_nodes, self.LS_GREY_DASHED);
                gplot(p_Adj, nodes_displaced, self.LS_BLACK_SOLID);
                hold off;

            end

        end

        % CALCULATE OPTIMUM ROW AND COLUMN NUMBER OF SUBPLOTS
        %
        % Inputs:
        %   p_N - overall number of plots
        function [rows, cols] = calcRowsCols(self, p_N)

            cols = ceil(sqrt(p_N));
            rows = ceil(p_N/cols);

        end

    end
end
