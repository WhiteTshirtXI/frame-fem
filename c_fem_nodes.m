classdef c_fem_nodes < handle
%C_FEM_NODES - Class containing node informations of the FEM system
% This class contains properties and functions for the handling of the
% FEM system nodes. Nodes can be added, defined, deleted and boundary
% conditions, nodal sources etc. can be applied.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_fem_nodes.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Constant properties :
%    NDOF - number of DOFs per node
%
% Properties :
%    x   - row vector of node x-coordinates
%    z   - row vector of node z-coordinates
%
%    bc  - row vector of nodal boundary conditions (fixed degrees of
%          freedom)
%    s   - row vector of nodal sources (forces, charge)
%
%    iln - last added node index
%
% Methods :
%    c_fem_nodes - constructor
%
%    addNode     - add node(s) and return indices
%    xz          - return node xz-coordinate vectors
%    n           - return number of nodes
%
%    fixNode     - fix nodal dof(s)
%    freeDOFs    - return free dofs mask vector 
%    idxFreeDOFs - return index vector of free dofs
%    setSource   - set nodal source(s)
%
%    distAl      - return distance and angle between two nodes
%
% Static Methods :
%    sysIdx      - return system dof indices for node(s)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-06-21 12:24 CEST
% Last Modified : 2012-08-10 11:24 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CONSTANT PROPERTIES %
    properties (Constant)

        % number of DOFs per node
        NDOF = 4;

    end


    % PRIVATE PROPERTIES %
    properties (SetAccess = private)

        % row vector of node x-coordinates
        x;
        % row vector of node z-coordinates
        z;

        % row vector of nodal boundary conditions (fixed degrees of
        % freedom)
        % 1 -> DOF is fixed
        % 0 -> DOF is free
        bc;

        % row vector of nodal sources (forces, moments, electric
        % potential)
        s;

        % last added node index
        iln;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_n - (maximum) number of nodes
        function s = c_fem_nodes( p_n )

            % pre allocate coordinate vectors
            s.x = sparse( 1, p_n );
            s.z = sparse( 1, p_n );

            % initialize boundary conditions and nodal sources as zero
            % matrices (each row corresponds to a node, each column
            % corresponds to a nodal degree of freedom)
            s.bc = sparse( s.NDOF, p_n );
            s.s  = sparse( s.NDOF, p_n );

            % last node index = 0
            s.iln = 0;

        end

        % ADD NODE(S) AND RETURN INDICES
        %
        % Inputs:
        %   p_x - node x-coordinate(s)
        %   p_z - node z-coordinate(s)
        function idx = addNode( s, p_x, p_z )

            % number of new nodes
            nNew = numel(p_x);

            % index list of newly added nodes
            idx = s.iln+1:s.iln+nNew;

            % index corrector for double node entries
            idx_corr = 0;

            % check if nodes already exist
            for i = 1:nNew

                dbl = find((s.x == p_x(i)) & (s.z == p_z(i)), 1);

                if (~isempty(dbl)) && (dbl <= s.iln)

                    % switch index
                    idx(i) = dbl;

                    % increment index corrector
                    idx_corr = idx_corr+1;

                else

                    % correct indices
                    idx(i) = idx(i)-idx_corr;

                end

            end

            % add node coordinates to coordinate vectors
            s.x(idx) = p_x(:).';
            s.z(idx) = p_z(:).';

            % update index of last added node
            s.iln = s.iln+nNew-idx_corr;

        end

        % RETURN NODE XZ-COORDINATE VECTORS
        %
        % Inputs:
        %   p_i - node index or index vector (OPTIONAL)
        function xz = xz( s, p_i )

            if nargin < 2
                p_i = 1:s.iln;
            end

            xz = [ s.x(p_i) ; s.z(p_i) ].';

        end

        % RETURN NUMBER OF NODES
        function n = n( s )

            n = numel(s.x);

        end

        % FIX NODAL DOF(s)
        %
        % Inputs:
        %   p_i - node index or index vector
        %   p_f - fixed dofs (1 -> dof fixed, 0 -> dof free)
        function s = fixNode( s, p_i, p_f )

            s.bc(:,p_i) = p_f;

        end

        % RETURN FREE DOFS MASK VECTOR
        function idx = freeDOFs( s )

            % all node DOFs index vector
            idx = ~s.bc(:);

        end

        % RETURN INDEX VECTOR OF FREE DOFS
        function idx = idxFreeDOFs( s )

            % all node DOFs index vector
            idx = (1:s.NDOF*s.n).';

            % use bc vector as mask for the index-vector
            idx = idx(s.freeDOFs());

        end

        % SET NODAL SOURCE(S)
        %
        % Inputs:
        %  p_i - node index or index vector
        %  p_s - source vector or source matrix (for multiple nodes)
        function s = setSource( s, p_i, p_s )

            s.s(:,p_i) = p_s;

        end

        % RETURN DISTANCE AND ANGLE BETWEEN TWO NODES
        %
        % Inputs:
        %   p_i - vector containing both node indices
        function [d,a] = distAl( s, p_i )

            xz12 = s.xz( p_i );

            % resulting vector
            v_xz12 = xz12(2,:) - xz12(1,:);

            % calculate distance and angle
            d = norm(v_xz12);
            a = -atan2(v_xz12(2), v_xz12(1));

        end

    end

    % STATIC METHODS
    methods (Static)
        
        % RETURN SYSTEM DOF INDICES FOR NODE(S)
        %
        % Inputs:
        %   p_i - node index or index vector
        function idx = sysIdx( p_i )

            NDOF = s.NDOF;
            
            % first: multiply the element node index vector by NDOF and
            %        substract (NDOF-1) to get the index vector in terms
            %        of global DOF indices
            idx_g = NDOF*p_i(:).'-(NDOF-1);
            
            % second: expand the node DOF index row vector idx_g to get
            %         a matrix with NDOF rows
            idx = repmat(idx_g, NDOF, 1);

            % third: expand the local DOF index vector [0;1;2] in the
            %        same way and sum up both matrices
            idx = idx + repmat([0;1;2], 1, numel(p_i));

            % fourth: the final dof index vector is generated by
            %         concatenating all of the matrix' columns
            idx = idx(:);

        end

        
    end


end
