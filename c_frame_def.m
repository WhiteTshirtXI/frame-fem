classdef c_frame_def < handle
%C_FRAME_DEF - Class for the definition of a 2D frame-structure
% This class is used to construct a twodimensional frame structure using
% beam-truss-like elements. These "parent"-elements however are
% subdivided into "child"-elements, which are used to construct a
% network of nodes and finite elements in the process of fem
% discretization.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_frame_def.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Properties :
%    IDXB_FIRSTNODE - beam vector index of first beam node
%    IDXB_LASTNODE  - beam vector index of last beam node
%    IDXB_RHOA      - beam vector index of rho*A
%    IDXB_EA        - beam vector index of E*A
%    IDXB_EI        - beam vector index of E*I
%    IDXB_NELEMENTS - beam vector index of number of beam elements
%
%
%    frame_nodes    - frame nodes array
%    frame_beams    - frame beams array
%
% Methods :
%    c_frame_def - constructor
%    addBeam     - add beam(s) to the frame structure
%    discretize  - discretize the frame structure using fem and return
%                  sys_fem-class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-25 11:05 CEST
% Last Modified : 2012-05-25 14:35 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)

        % beams array indices
        IDXB_FIRSTNODE = 1;
        IDXB_LASTNODE  = 2;
        IDXB_RHOA      = 3;
        IDXB_EA        = 4;
        IDXB_EI        = 5;
        IDXB_NELEMENTS = 6;

        % frame nodes array ([x1 z1; x2 z2; ...])
        frame_nodes;

        % frame beams array ( [node1, node2, rhoA, EA, EI, nElements] )
        frame_beams = [];

    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_nodes - x- and z-coordinate vector of the frame nodes
        function self = c_frame_def( p_nodes )

            self.frame_nodes = p_nodes;

        end

        % ADD BEAM(S) TO THE FRAME STRUCTURE
        %
        % Inputs:
        %   p_beams - specification vectors for the beams. each row
        %             contains a beam spec-vector with the following
        %             data:
        %             [ startNodeID endNodeID rhoA EA EI nElements ]
        function self = addBeam( self, p_beams )

            self.frame_beams = [self.frame_beams ; p_beams];

        end

        % DISCRETIZE THE FRAME STRUCTURE USING FEM AND RETURN
        % SYS_FEM-CLASS
        function sys_fem = discretize( self )

            % calculate number of system nodes
            N = size(self.frame_nodes,1) + ...
                sum(self.frame_beams(:,self.IDXB_NELEMENTS)-1);

            % allocate system nodes-vector
            sys_nodes = zeros(N, 2);

            % this vector contains the information of frame nodes, which
            % have already been added to the system nodes-vector and
            % what index they have
            % an index of 0 means, the frame node hasn't been assigned
            % already
            frame_nodes_index = zeros(size(self.frame_nodes,1),1);

            % node counter, starting with 1
            i_node = 1;

            % calculate nodes coordinates for each beam
            for beam = self.frame_beams'

                % beam first and last frame node index
                beamFirstNode = beam(self.IDXB_FIRSTNODE);
                beamLastNode = beam(self.IDXB_LASTNODE);

                % beam number of elements
                beamNEl = beam(self.IDXB_NELEMENTS);
                
                % start and end node coordinates
                xzStart = self.frame_nodes(beamFirstNode,:);
                xzEnd   = self.frame_nodes(beamLastNode,:);
                
                % check if the first frame node already hasn't been
                % assigned to the node indices of the discretized system
                if frame_nodes_index(beamFirstNode) == 0
                    % add start node coordinates to system nodes vector
                    sys_nodes(i_node,:) = xzStart;
                    % mark first beam node as 'assigned'
                    frame_nodes_index(beamFirstNode) = i_node;
                    % increment node counter by one
                    i_node = i_node + 1;
                end

                % beam vector between start and end node
                xzBeam = xzEnd - xzStart;

                % get beam angle and beam length
                lBeam = norm(xzBeam);
                alphaBeam = -atan2(xzBeam(2),xzBeam(1));

                % run through sub-nodes for the current beam
                for i_subNode = 1:(beamNEl-1)

                    % calculate sub-node-coordinates
                    xzSubNode = xzStart + i_subNode*xzBeam/beamNEl;

                    % add sub-node to system nodes vector
                    sys_nodes(i_node,:) = xzSubNode;

                    % increment node counter by one
                    i_node = i_node + 1;

                end
                
                % check if the last frame node already hasn't been
                % assigned to the node indices of the discretized system
                if frame_nodes_index(beamLastNode) == 0
                    % add start node coordinates to system nodes vector
                    sys_nodes(i_node,:) = xzEnd;
                    % mark first beam node as 'assigned'
                    frame_nodes_index(beamLastNode) = i_node;
                    % increment node counter by one
                    i_node = i_node + 1;
                end

            end

            % initialize fem-system class using the node list
            sys_fem = c_sys_fem(sys_nodes);

        end
    end

end
