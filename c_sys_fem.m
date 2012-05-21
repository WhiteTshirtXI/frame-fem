classdef c_sys_fem < handle
%C_SYS_FEM - Class representing a FEM system
% This class contains properties and methods for the discretization of
% a system of specified nodes and elements via FEM. The resulting system
% matrices can easily be accessed by the access-methods.
%
% ! currently only beam-truss-like elements are supported !
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_sys_fem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Properties :
%    NDOF      - number of DOFs per node
%    IDOF_U    - DOF index for x-translation
%    IDOF_W    - DOF index for z-translation
%    IDOF_B    - DOF index for rotation around y-axis
%
%    IDOF_VEC  - node DOF index vector
%
%    N         - number of nodes
%    nodes     - node list
%
%    MSys      - system mass matrix
%    KSys      - systen stiffness matrix
%
%    eigVal    - matrix of eigenvalues of the system
%    eigVec    - matrix of eigenvectors of the system
%    eigRecalc - recalculation of eigenvalues required
%
% Methods :
%    c_sys_fem   - constructor
%    sysDOF      - return overall system degrees of freedom
%    add_element - add beam element
%    eigCalc     - calculate eigenvalues and eigenvectors
%    getModes    - return eigenmodes
%    getAllModes - return all eigenmodes
%    removeModes - remove eigenmodes with an eigenfrequency below a
%                  certain tolerance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-18 12:50 CEST
% Last Modified : 2012-05-21 13:09 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)
        % number of DOFs per node
        % (equals 3, because, currently, only beam-truss-like elements
        % are supported)
        NDOF = 3;

        % DOF indices (currently unused)
        %  u    -> 1
        %  w    -> 2
        %  beta -> 3
        IDOF_U = 1;
        IDOF_W = 2;
        IDOF_B = 3;

        % node DOF index vector
        IDOF_VEC;

        % overall number of nodes
        N = 0;

        % node list (format: [x1 z1 ; x2 z2 ; ...])
        nodes;

        % system mass and stiffness matrices
        MSys;
        KSys;

        % eigenvalues and -vectors
        eigVal;
        eigVec;

        % recalculation of eigenvalues required
        eigRecalc = true;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        %
        % Inputs:
        %   p_nodes - x- and z-coordinate vector of the system nodes
        function self = c_sys_fem( p_nodes )

            self.nodes = p_nodes;
            % the number of system nodes equals the number of rows of
            % the array p_nodes
            self.N     = size(p_nodes, 1);

            % overall degrees of freedom
            sysDOF = self.sysDOF();

            % create node DOF index vector
            self.IDOF_VEC = [1:self.NDOF];

            % allocate memory for the system matrices
            self.MSys  = sparse(sysDOF,sysDOF);
            self.KSys  = sparse(sysDOF,sysDOF);

            % allocte memory for the eigenvalue matrices
            % self.eigVal = spalloc(sysDOF,sysDOF,sysDOF);
            self.eigVal = zeros(sysDOF);
            self.eigVec = zeros(sysDOF);

        end

        % RETURN OVERALL SYSTEM DEGREES OF FREEDOM
        function sysDOF = sysDOF(self)

            sysDOF = self.N * self.NDOF;

        end

        % ADD BEAM ELEMENT
        %
        % Inputs:
        %   p_idx  - index vector matrix for local/global
        %          - node-assigmnent [ idx1_el1 idx2_el1 ;
        %                              idx1_el2 idx2_el2 ; ... etc.
        %   p_rhoA - element line density
        %   p_EA   - element tensile rigidity
        %   p_EI   - element bending rigidity
        function self = add_element(self, p_idx, p_rhoA, p_EA, p_EI)

            % run through the rows of p_idx
            for idx = p_idx'

                % get x- and z-coordinates of both element nodes using
                % the index vector
                xz12 = self.nodes(idx',:);

                % resulting vector between both nodes
                dxz12 = xz12(2,:) - xz12(1,:);

                % calculate element length and angle
                le    = norm(dxz12);
                alpha = -atan2(dxz12(2), dxz12(1));

                % calculate element matrices
                [ Me, Ke ] = matrices_beam( le, alpha, ...
                                            [p_rhoA p_EA p_EI] );

                % calculate DOF index vector for matrix summation to
                % system matrix
                %
                % first: multiply the element node index vector by NDOF
                %        and substract (NDOF-1) to get the index vector
                %        in terms of global DOF indices
                g_idx_dof = self.NDOF*idx'-(self.NDOF-1);
                
                % second: expand the node DOF index column vector
                %         p_idx_dof to get a matrix with NDOF columns
                %         and two (number of element nodes) rows
                m_idx_dof = repmat(g_idx_dof', 1, self.NDOF);

                % third: expand the local DOF index vector IDOF_VEC-1 in
                %        the same way and sum up both matrices
                m_idx_dof = m_idx_dof + repmat(self.IDOF_VEC-1, 2, 1);

                % fourth: the final dof index vector is generated by
                %         concatenating all of the matrix' columns
                idx_dof = m_idx_dof(:);

                % using these index vectors, add the element matrices to
                % the system matrices
                self.MSys(idx_dof,idx_dof) = self.MSys(idx_dof, ...
                                                       idx_dof) + Me;    
                self.KSys(idx_dof,idx_dof) = self.KSys(idx_dof, ...
                                                       idx_dof) + Ke;    

            end

            % the system matrices have changed, so a recalculation of the
            % eigenvalues is required
            self.eigRecalc = true;

        end 

        % CALCULATE EIGENVALUES
        function self = eigCalc(self)

            % check if eigenvalue recalculation is neccessary
            if(self.eigRecalc)

                % WORKAROUND : using full() to convert the sparse system
                %              matrices to full matrices
                [self.eigVec,self.eigVal] = eig( full(self.KSys), ...
                                                -full(self.MSys));

                self.eigRecalc = false;

            end

        end

        % RETURN EIGENMODES
        %
        % Inputs:
        %   p_nModes - number of modes to be returned
        function [lambda,V] = getModes(self, p_nModes)

            % if neccessary, recalc the eigenmodes
            self.eigCalc();

            % clean up eigenmodes
            [lambda_c,V_c,nDel] = self.removeModes();

            % the eigenvalues are fetched in the order from lowest value
            % to highest value. the lowest eigenvalue is located in the
            % last cell of the eigenvalue-matrix, so the index of the
            % highest eigenvalue needs to be calculated according to
            % p_nModes.
            idx_maxLambda = max(self.sysDOF()-nDel-p_nModes+1, 1);

            % eigenvalues vector
            lambda = flipud(diag(lambda_c(idx_maxLambda:end, ...
                                          idx_maxLambda:end)));

            % corresponding eigenvectors
            V = V_c(:,idx_maxLambda:end);

        end

        % RETURN ALL EIGENMODES
        function [lambda,V] = getAllModes(self)

            [lambda,V] = self.getModes(self.sysDOF());

        end

        % REMOVE EIGENMODES WITH AN EIGENFREQUENCY BELOW A CERTAIN
        % TOLERANCE
        function [lambda,V,nDel] = removeModes(self)

            % minimum frequency
            OMEGA_MIN = 1*2*pi;

            % calculate eigenfrequencies
            freq = imag(sqrt(diag(self.eigVal)));

            % create mask vector by checking if the frequencies are
            % above the speciefied frequency threshold
            mask = freq > OMEGA_MIN;

            % remove the affected rows and columns using the mask vector
            lambda = self.eigVal(mask,mask);
            % in the eigenvector matrix only the affected colors have to
            % be removed
            V = self.eigVec(:,mask);

            % number of deleted elements
            nDel = nnz(~mask);

        end

    end

end
