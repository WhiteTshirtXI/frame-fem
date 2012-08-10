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
%    frame_nodes     - frame nodes array
%    frame_beams     - frame beams array
%    frame_nodes_idx - frame nodes index vector for the discretized
%                      system
%    frame_nodes_bc  - frame nodal boundary conditions array
%    frame_inf_bc    - frame nodal boundary conditions array for
%                      infinite elements
%    frame_sh        - nodal harmonic source amplitudes vector
%
% Methods :
%    c_frame_def       - constructor
%
%    addBeam           - add beam(s) to the frame structure
%    discretize        - discretize the frame structure using fem and
%                        return sys_fem-class
%
%    nodeBC_clamped    - apply clamped boundary condition to frame nodes
%    nodeBC_spprted    - apply simply supported boundary condition to
%                        frame nodes
%    nodeBC_infinite   - apply infinite element boundary condition to
%                        frame nodes (WANG & LAI 1999)
%
%    addHarmonicSource - add harmonic source to nodes
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
% Last Modified : 2012-08-09 15:20 CEST
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

        % frame nodes index vector for the discretized system
        frame_nodes_idx;

        % frame nodal boundary conditions array
        frame_nodes_bc;

        % frame infinite element bcs
        % contains stiffness and damping coefficients for each node and
        % nodal DOF (if not zero)
        frame_inf_bc;

        % nodal harmonic source amplitudes vector
        % ([fx1 fz1 mxz1 sigma1; fx2 fz2 mxz2 sigma2; ...])
        frame_sh;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_nodes - x- and z-coordinate vector of the frame nodes
        function s = c_frame_def( p_nodes )

            s.frame_nodes = p_nodes;
    
            % initialize nodes index vector as empty
            s.frame_nodes_idx = zeros(size(s.frame_nodes,1),1);

            % initialize boundary conditions vectors as empty
            s.frame_nodes_bc = zeros(size(p_nodes,1),4);
            s.frame_inf_bc = zeros(size(p_nodes,1),8);

            % initialize nodal sources as empty
            s.frame_sh = zeros(size(p_nodes,1),4);

        end

        % ADD BEAM(S) TO THE FRAME STRUCTURE
        %
        % Inputs:
        %   p_beams - specification vectors for the beams. each row
        %             contains a beam spec-vector with the following
        %             data:
        %             [ startNodeID endNodeID rhoA EA EI elDensity ]
        function s = addBeam( s, p_beams )

            % calculate number of elements on each beam via the
            % specified element density

            % beam node indices
            bIdx = p_beams(:,[s.IDXB_FIRSTNODE s.IDXB_LASTNODE]);

            % beam lengths in x- and z-direction as a difference between
            % the last and first node x- and z-coordinate
            bLx = s.frame_nodes(bIdx(:,2),1)- ...
                  s.frame_nodes(bIdx(:,1),1);
            bLz = s.frame_nodes(bIdx(:,2),2)- ...
                  s.frame_nodes(bIdx(:,1),2);

            % number of elements = beam length * element density
            p_beams(:,s.IDXB_NELEMENTS) = ceil( ...
                   sqrt(bLx.^2+bLz.^2).*p_beams(:,s.IDXB_NELEMENTS));


            s.frame_beams = [s.frame_beams ; p_beams];

        end

        % DISCRETIZE THE FRAME STRUCTURE USING FEM AND RETURN
        % SYS_FEM-CLASS
        function sys_fem = discretize( s )

            % re-initialize frame nodes index vector
            s.frame_nodes_idx = zeros(size(s.frame_nodes,1),1);

            % calculate number of system nodes
            N = size(s.frame_nodes,1) + ...
                sum(s.frame_beams(:,s.IDXB_NELEMENTS)-1);

            % calculate number of system elements
            NEl = sum(s.frame_beams(:,s.IDXB_NELEMENTS));

            % initialize fem-system class
            sys_fem = c_sys_fem(N, NEl);

            % for each frame beam add nodes and elements to the
            % fem-system class
            for beam = s.frame_beams.'

                % beam first and last frame node index
                beamFirstNode = beam(s.IDXB_FIRSTNODE);
                beamLastNode = beam(s.IDXB_LASTNODE);

                % beam number of elements
                beamNEl = beam(s.IDXB_NELEMENTS);

                % beam geometrical and material properties
                beamRhoA = beam(s.IDXB_RHOA);
                beamEA = beam(s.IDXB_EA);
                beamEI = beam(s.IDXB_EI);
                
                % start and end node coordinates
                xStart = s.frame_nodes(beamFirstNode,1);
                xEnd   = s.frame_nodes(beamLastNode,1);
                zStart = s.frame_nodes(beamFirstNode,2);
                zEnd   = s.frame_nodes(beamLastNode,2);

                % length and orientation angle
                l = sqrt((xEnd-xStart)^2+(zEnd-zStart)^2);
                a = -atan2(zEnd-zStart,xEnd-xStart);

                % beam node coordinates
                xNodes = linspace(xStart,xEnd,beamNEl+1);
                zNodes = linspace(zStart,zEnd,beamNEl+1);

                % add nodes to fem-system class and get indices
                nodeIdx = sys_fem.addNodes(xNodes, zNodes);

                % insert beam elements with constant properties between
                % these nodes
                sys_fem.addElBeamConst(nodeIdx, l/beamNEl, a,       ...
                                       beamRhoA, beamEA, beamEI);

                % update frame nodes indices
                s.frame_nodes_idx([beamFirstNode beamLastNode]) = ...
                                                       nodeIdx([1 end]);

            end % beam loop

            % apply nodal boundary conditions
            sys_fem.addNodeBC(s.frame_nodes_idx, ...
                              s.frame_nodes_bc);

            % apply nodal boundary conditions for infinite elements
            sys_fem.addNodeBC_inf(s.frame_nodes_idx, ...
                                  s.frame_inf_bc);

            % add harmonic source amplitudes to system
            sys_fem.addNodeSources(s.frame_nodes_idx, s.frame_sh);

        end

        % APPLY CLAMPED BOUNDARY CONDITION TO FRAME NODES
        %
        % Inputs:
        %   p_nodes - node indices of the nodes which need to be clamped
        function s = nodeBC_clamped( s, p_nodes )

            for node = p_nodes
                s.frame_nodes_bc(node,:) = [1 1 1 0];
            end

        end

        % APPLY SIMPLY SUPPORTED BOUNDARY CONDITION TO FRAME NODES
        %
        % Inputs:
        %   p_nodes - node indices of the nodes which need to be simply
        %             supported
        function s = nodeBC_spprted( s, p_nodes )

            for node = p_nodes
                s.frame_nodes_bc(node,:) = [1 1 0 0];
            end

        end

        % APPLY INFINITE ELEMENT BOUNDARY CONDITION TO FRAME NODES
        % (WAND & LAI 1999)
        %
        % Inputs:
        %   p_nodes - node indices of the nodes where the infinite
        %             elements need to be attached
        %   p_omega - angular frequency of excitation
        function s = nodeBC_infinite( s, p_nodes, p_omega )
            
            % clear boundary condition vector
            s.frame_inf_bc = zeros(size(s.frame_nodes,1),8);

            for node = p_nodes

                % TODO: this does not account for nodes, where multiple
                % beams are connected to!
                % TODO: change for arbitrarily oriented beams!
                
                % find firstmost beam ID which is connected to the
                % specified node
                [beamID,~] = find(s.frame_beams(:,[s.IDXB_FIRSTNODE  ...
                                           s.IDXB_LASTNODE]) == node,1);

                % calculate wave numbers for the beam
                beamRhoA = s.frame_beams(beamID,s.IDXB_RHOA);
                beamEA   = s.frame_beams(beamID,s.IDXB_EA);
                beamEI   = s.frame_beams(beamID,s.IDXB_EI);

                % longitudinal wave number
                kL = p_omega*sqrt(beamRhoA/beamEA);
                % bending wave number
                kB = sqrt(p_omega)*(beamRhoA/beamEI)^0.25;

                % calculate spring stiffness and damping coefficient for
                % longitudinal case
                bcKL = 0;
                bcDL = beamEA*kL/p_omega;
                % calculate spring stiffness and damping coefficient for
                % bending case
                bcKB = 0.5*beamEI*kB^3;
                bcDB = bcKB/p_omega;

                % add to boundary condition vector
                s.frame_inf_bc(node,:) = [ bcKL bcDL bcKB bcDB 0 0 0 0];

            end

        end

        % ADD HARMONIC SOURCE TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force is applied
        %   p_S      - source amplitudes magnitude,
        %              [ Fx ; Fz ; My ; sigma ]
        function s = addHarmonicSource( s, p_nodes, p_S )

            for node = p_nodes
                s.frame_sh(node,:) = p_S.';
            end

        end

        % ADD HARMONIC FORCE TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force is applied
        %   p_F      - force amplitudes magnitude, [ Fx ; Fz ; My ]
        function s = addHarmonicForce( s, p_nodes, p_F )

            s.addHarmonicSource(p_nodes, [p_F ; 0]);

        end

        % ADD HARMONIC FORCE IN X-DIRECTION TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force in x-direction is applied
        %   p_F      - force amplitude magnitude
        function s = addHarmonicForceX( s, p_nodes, p_F )

            s.addHarmonicForce(p_nodes, [p_F ; 0 ; 0]);

        end

        % ADD HARMONIC FORCE IN Z-DIRECTION TO NODES
        %
        % Inputs:
        %   p_nodes  - node indices of the nodes where an external
        %              harmonic force in z-direction is applied
        %   p_F      - force amplitude magnitude
        function s = addHarmonicForceZ( s, p_nodes, p_F )

            s.addHarmonicForce(p_nodes, [0 ; p_F ; 0]);

        end

    end

end
