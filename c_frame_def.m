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
%    frame_nodes_bc - frame nodal boundary conditions array
%    frame_inf_bc   - frame nodal boundary conditions array for infinite 
%                     elements
%    frame_fh       - nodal harmonic force amplitudes vector
%
% Methods :
%    c_frame_def       - constructor
%
%    addBeam           - add beam(s) to the frame structure
%    discretize        - discretize the frame structure using fem and
%                        return sys_fem-class
%
%    nodeBC_clamped    - apply clamped boundary condition to frame nodes
%    nodeBC_jointed    - apply jointed boundary condition to frame nodes
%    nodeBC_infinite   - apply infinite element boundary condition to
%                        frame nodes (WANG & LAI 1999)
%
%    addHarmonicForce  - add harmonic force to nodes
%    addHarmonicForceX - add harmonic force in x-direction to nodes
%    addHarmonicForceZ - add harmonic force in z-direction to nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-25 11:05 CEST
% Last Modified : 2012-05-31 17:24 CEST
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
        frame_beams;

        % frame nodal boundary conditions array
        frame_nodes_bc;

        % frame infinite element bcs
        % contains [ K_bc D_bc ] for each node (if not zero)
        frame_inf_bc;

        % nodal harmonic force amplitudes vector
        % ([fx1 fz1 mxz1; fx2 fz2 mxz2; ...])
        frame_fh;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_nodes - x- and z-coordinate vector of the frame nodes
        function self = c_frame_def( p_nodes )

            self.frame_nodes = p_nodes;

            % initialize boundary conditions vectors as empty
            self.frame_nodes_bc = zeros(size(p_nodes,1),3);
            self.frame_inf_bc = zeros(size(p_nodes,1),2);

            % initialize nodal forces and moments as empty
            self.frame_fh = zeros(size(p_nodes,1),3);

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

            % calculate number of system elements
            NEl = sum(self.frame_beams(:,self.IDXB_NELEMENTS));

            % allocate system nodes-vector
            sys_nodes = zeros(N, 2);

            % allocate system elements-cell array
            % each row contains: { [elNode1 elNode2], [rhoA EA EI] }
            sys_elements = cell(NEl, 2);

            % this vector contains the information of frame nodes, which
            % have already been added to the system nodes-vector and
            % what index they have
            % an index of 0 means, the frame node hasn't been assigned
            % already
            frame_nodes_index = zeros(size(self.frame_nodes,1),1);

            % node counter, starting with 1
            i_node = 1;

            % element counter, starting with 1
            i_element = 1;

            % calculate nodes coordinates for each beam
            for beam = self.frame_beams'

                % beam first and last frame node index
                beamFirstNode = beam(self.IDXB_FIRSTNODE);
                beamLastNode = beam(self.IDXB_LASTNODE);

                % beam number of elements
                beamNEl = beam(self.IDXB_NELEMENTS);

                % beam geometrical and material properties
                beamRhoA = beam(self.IDXB_RHOA);
                beamEA = beam(self.IDXB_EA);
                beamEI = beam(self.IDXB_EI);
                
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

                % first element node index
                elNode1 = frame_nodes_index(beamFirstNode);

                % run through sub-nodes for the current beam
                for i_subNode = 1:(beamNEl-1)

                    % calculate sub-node-coordinates
                    xzSubNode = xzStart + i_subNode*xzBeam/beamNEl;

                    % add sub-node to system nodes vector
                    sys_nodes(i_node,:) = xzSubNode;

                    % second element node index
                    elNode2 = i_node;

                    % increment node counter by one
                    i_node = i_node + 1;

                    % add element data to system elements-cell vector
                    sys_elements(i_element,:) = {[elNode1 elNode2] ...
                                                 [beamRhoA beamEA ...
                                                  beamEI] };

                    % increment element counter by one
                    i_element = i_element + 1;

                    % switch first and second element node index
                    elNode1 = elNode2;

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

                % second element node index
                elNode2 = frame_nodes_index(beamLastNode);

                % add element data to system elements-cell vector
                sys_elements(i_element,:) = {[elNode1 elNode2] ...
                                             [beamRhoA beamEA beamEI] };

                % increment element counter by one
                i_element = i_element + 1;

            end

            % initialize fem-system class using the node list
            sys_fem = c_sys_fem(sys_nodes);

            % apply nodal boundary conditions
            sys_fem.addNodeBC(frame_nodes_index, self.frame_nodes_bc);

            % apply nodal boundary conditions for infinite elements
            sys_fem.addNodeBC_inf(frame_nodes_index, ...
                                  self.frame_inf_bc);

            % add harmonic force amplitudes to system
            sys_fem.addNodeForces(frame_nodes_index, self.frame_fh);

            % buld force vector
            sys_fem.buildFSys();

            % add all beam elements to the fem system
            sys_fem.add_elements(sys_elements);

        end

        % APPLY CLAMPED BOUNDARY CONDITION TO FRAME NODES
        %
        % Inputs:
        %   p_nodes - node indices of the nodes which need to be clamped
        function self = nodeBC_clamped( self, p_nodes )

            for node = p_nodes
                self.frame_nodes_bc(node,:) = [1 1 1];
            end

        end

        % APPLY JOINTED BOUNDARY CONDITION TO FRAME NODES
        %
        % Inputs:
        %   p_nodes - node indices of the nodes which need to be joint
        function self = nodeBC_jointed( self, p_nodes )

            for node = p_nodes
                self.frame_nodes_bc(node,:) = [1 1 0];
            end

        end

        % APPLY INFINITE ELEMENT BOUNDARY CONDITION TO FRAME NODES
        % (WAND & LAI 1999)
        %
        % Inputs:
        %   p_nodes - node indices of the nodes where the infinite
        %             elements need to be attached
        %   p_omega - angular frequency of excitation
        function self = nodeBC_infinite( self, p_nodes, p_omega )
            
            % clear boundary condition vector
            self.frame_inf_bc = zeros(size(self.frame_nodes,1),2);

            for node = p_nodes

                % TODO: this does not account for nodes, where multiple
                % beams are connected to!
                % TODO: change for arbitrarily oriented beams!
                
                % find firstmost beam ID which is connected to the
                % specified node
                [beamID,beamIDcol] = find(                           ...
                   self.frame_beams(:,[self.IDXB_FIRSTNODE          ...
                                       self.IDXB_LASTNODE]) == node,1);

                % calculate bending wave number for the beam
                beamRhoA = self.frame_beams(beamID,self.IDXB_RHOA);
                beamEI   = self.frame_beams(beamID,self.IDXB_EI);

                kB = sqrt(p_omega)*(beamRhoA/beamEI)^0.25;

                % calculate spring stiffness and damping coefficient
                bcK = 0.5*beamEI*kB^3;
                bcD = bcK/p_omega;

                % add to boundary condition vector
                self.frame_inf_bc(node,:) = [ bcK bcD ];

            end

        end

        % ADD HARMONIC FORCE TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force is applied
        %   p_F      - force amplitude magnitude
        %   p_Falpha - force direction (    0 -> pos. x-Direction
        %                               +pi/2 -> neg. z-Direction
        %                                  pi -> neg. x-Direction
        %                               -pi/2 -> pos. z-Direction )
        function self = addHarmonicForce( self, p_nodes, p_F, p_Falpha )

            fx =  p_F*cos(p_Falpha);
            fz = -p_F*sin(p_Falpha);

            for node = p_nodes
                self.frame_fh(node,1:2) = [fx fz];
            end

        end

        % ADD HARMONIC FORCE IN X-DIRECTION TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force in x-direction is applied
        %   p_F      - force amplitude magnitude
        function self = addHarmonicForceX( self, p_nodes, p_F )

            self.addHarmonicForce(p_nodes, p_F, 0);

        end

        % ADD HARMONIC FORCE IN Z-DIRECTION TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force in z-direction is applied
        %   p_F      - force amplitude magnitude
        function self = addHarmonicForceZ( self, p_nodes, p_F )

            self.addHarmonicForce(p_nodes, p_F, -pi/2);

        end

    end

end
