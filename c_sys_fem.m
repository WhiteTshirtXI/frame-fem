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
%    NDOF       - number of DOFs per node
%    IDOF_U     - DOF index for x-translation
%    IDOF_W     - DOF index for z-translation
%    IDOF_B     - DOF index for rotation around y-axis
%
%    IDOF_VEC   - node DOF index vector
%
%    N          - number of nodes
%    nodes      - node list
%    nAdj       - nodes adjacancy matrix
%
%    bc         - boundary conditions for each node and nodal DOF
%    bc_inf     - infinite boundary conditions for each node
%    fh         - harmonic force amplitudes for each node
%
%    MSys       - system mass matrix
%    KSys       - system stiffness matrix
%    fSys       - system nodal force vector
%
%    eigVal     - matrix of eigenvalues of the system
%    eigVec     - matrix of eigenvectors of the system
%    eigRecalc  - recalculation of eigenvalues required
%
%    plot_nodes - nodes plotting class
%
% Methods :
%    c_sys_fem      - constructor
%    sysDOF         - return overall system degrees of freedom
%    add_element    - add beam element
%    add_elements   - add multiple beam elements
%
%    eigCalc        - calculate eigenvalues and eigenvectors
%    eigOmega       - calculate the angular eigenfrequencies
%    eigF           - calculate the eigenfrequencies
%
%    getModes       - return eigenmodes
%    getAllModes    - return all eigenmodes
%    removeModes    - remove eigenmodes with an eigenfrequency below a
%                     certain tolerance
%    plotModes      - plot eigenmodes
%    addNodeBC      - add nodal boundary condition
%    addNodeBC_inf  - add nodal boundary condition for infinite elements
%    nodeBC_clamped - clamped nodal BC (all three DOFs fixed)
%    nodeBC_jointed - jointed nodal BC (rotational DOF free)
%
%    getFixedDOFs   - return vector of fixed dofs
%
%    addNodeForces  - add nodal forces and moments
%    buildFSys      - build system force vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-18 12:50 CEST
% Last Modified : 2012-05-30 15:10 CEST
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

        % nodes adjacancy matrix
        nAdj;

        % boundary conditions for each node and nodal DOF
        bc;

        % infinite boundary conditions for each node
        bc_inf;

        % harmonic force amplitudes for each node
        fh;

        % system mass and stiffness matrices
        MSys;
        KSys;

        % system nodal force vector
        fSys;

        % eigenvalues and -vectors
        eigVal;
        eigVec;

        % recalculation of eigenvalues required
        eigRecalc = true;

        % nodes plotting class
        plot_nodes;

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

            % allocate memory for the adjacancy matrix
            self.nAdj = sparse(self.N,self.N);

            % initialize all boundary conditions to zero (non-fixed)
            self.bc = sparse(self.N,self.NDOF);

            % initialize infinite bcs to zero
            self.bc_inf = sparse(self.N,2);

            % initialize force amplitudes to zero
            self.fh = sparse(self.N,self.NDOF);

            % overall degrees of freedom
            sysDOF = self.sysDOF();

            % create node DOF index vector
            self.IDOF_VEC = [1:self.NDOF];

            % allocate memory for the system matrices and vectors
            self.MSys  = sparse(sysDOF,sysDOF);
            self.KSys  = sparse(sysDOF,sysDOF);
            self.fSys  = sparse(sysDOF,sysDOF);

            % allocate memory for the eigenvalue matrices
            % self.eigVal = spalloc(sysDOF,sysDOF,sysDOF);
            self.eigVal = zeros(sysDOF);
            self.eigVec = zeros(sysDOF);

            % create plot nodes class
            self.plot_nodes = c_plot_nodes(p_nodes);

        end

        % RETURN OVERALL SYSTEM DEGREES OF FREEDOM
        function sysDOF = sysDOF(self)

            sysDOF = self.N * self.NDOF;

        end

        % ADD BEAM ELEMENT
        %
        % Inputs:
        %   p_idx   - node indices [idx1 idx2]
        %   p_bProp - element properties [rhoA EA EI]
        function self = add_element(self, p_idx, p_bProp)

            % get x- and z-coordinates of both element nodes using the
            % index vector
            xz12 = self.nodes(p_idx,:);

            % resulting vector between both nodes
            dxz12 = xz12(2,:) - xz12(1,:);

            % calculate element length and angle
            le    = norm(dxz12);
            alpha = -atan2(dxz12(2), dxz12(1));

            % calculate element matrices
            [ Me, Ke ] = matrices_beam( le, alpha, p_bProp );

            % calculate DOF index vector for matrix summation to system
            % matrix
            %
            % first: multiply the element node index vector by NDOF and
            %        substract (NDOF-1) to get the index vector in terms
            %        of global DOF indices
            g_idx_dof = self.NDOF*p_idx'-(self.NDOF-1);
            
            % second: expand the node DOF index column vector p_idx_dof
            %         to get a matrix with NDOF columns and two (number
            %         of element nodes) rows
            m_idx_dof = repmat(g_idx_dof, 1, self.NDOF);

            % third: expand the local DOF index vector IDOF_VEC-1 in the
            %        same way and sum up both matrices
            m_idx_dof = m_idx_dof + repmat(self.IDOF_VEC-1, 2, 1);

            % fourth: the final dof index vector is generated by
            %         concatenating all of the matrix' columns
            idx_dof = m_idx_dof(:);

            % using these index vectors, add the element matrices to the
            % system matrices
            self.MSys(idx_dof,idx_dof) = self.MSys(idx_dof,idx_dof) + Me;
            self.KSys(idx_dof,idx_dof) = self.KSys(idx_dof,idx_dof) + Ke;


            % the system matrices have changed, so a recalculation of
            % the eigenvalues is required
            self.eigRecalc = true;

            % add entries to adjacancy matrix according to the node
            % indices
            self.nAdj(p_idx,p_idx) = ~eye(2);

        end 

        % ADD MULTIPLE BEAM ELEMENTS
        %
        % Inputs:
        %   p_beams - cell-array containing all the data necessary for
        %             the addition of every beam element to the system
        %             { [idx1 idx2] [rhoA EA EI] }
        function self = add_elements(self, p_beams)
            
            for beam = p_beams'

                self.add_element(beam{1},beam{2});

            end

        end

        % CALCULATE EIGENVALUES
        function self = eigCalc(self)

            % check if eigenvalue recalculation is neccessary
            if(self.eigRecalc)

                sysDOF = self.sysDOF();

                % mask vector for fixed DOF representation
                bcMask = ~self.getFixedDOFs();

                % number of non-fixed DOFs
                bcNFree = nnz(bcMask);

                % initialize eigenvector and eigenvalue matrices as zero
                self.eigVec = zeros(sysDOF,bcNFree);
                self.eigVal = zeros(bcNFree);

                % WORKAROUND : using full() to convert the sparse system
                %              matrices to full matrices
                [self.eigVec(bcMask,:),self.eigVal] =                ...
                                eig( full(self.KSys(bcMask,bcMask)), ...
                                    -full(self.MSys(bcMask,bcMask)));

                self.eigRecalc = false;

            end

        end

        % CALCULATE ANGULAR EIGENFREQUENCIES
        %
        % Inputs:
        %   p_nModes - number of eigenfrequencies to be returned (may be
        %              ommited to return all eigenfrequencies!)
        function omega = eigOmega(self, p_nModes)

            % check if the parameter p_nModes has been specified
            if exist('p_nModes')
                [lambda,V] = self.getModes(p_nModes);
            else
                [lambda,V] = self.getAllModes();
            end

            % return vector of angular eigenfrequencies
            omega = imag(sqrt(lambda));

        end

        % CALCULATE EIGENFREQUENCIES
        %
        % Inputs:
        %   p_nModes - number of eigenfrequencies to be returned (may be
        %              ommited to return all eigenfrequencies!)
        function f = eigF(self, p_nModes)

            % check if the parameter p_nModes has been specified
            if exist('p_nModes')
                f = self.eigOmega(p_nModes)./(2*pi);
            else
                f = self.eigOmega()./(2*pi);
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
            idx_maxLambda = max(size(self.eigVal,1)-nDel-p_nModes+1,1);

            % eigenvalues vector
            lambda = flipud(diag(lambda_c(idx_maxLambda:end, ...
                                          idx_maxLambda:end)));

            % corresponding eigenvectors
            V = fliplr(V_c(:,idx_maxLambda:end));

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

        % PLOT EIGENMODES
        %
        % Inputs:
        %   p_nModes - number of modes to be plotted
        function self = plotModes(self, p_nModes)

            % if neccessary, recalc the eigenmodes
            self.eigCalc();

            % get modes
            [lambda, V] = self.getModes(p_nModes);

            % number of acquired modes equals number of columns in V
            nModes = size(V, 2);

            % allocate nodal displacement matrix with Nx(2*nModes)
            % cells
            nodeDisplacements = zeros(self.N, 2*nModes);

            % assemble nodal displacement matrix
            % first: insert all x-displacements in the eigenvector to
            %        every odd column of the nodal displacement matrix
            nodeDisplacements(:,1:2:end) = V(1:self.NDOF:end,:);

            % second: insert all z-displacements in the eigenvector to
            %         every even column of the nodal displacement matrix
            nodeDisplacements(:,2:2:end) = V(2:self.NDOF:end,:);

            % create titles for each subplot
            
            % calculate frequency
            f = imag(sqrt(lambda))/(2*pi);
            % mode numbers
            mNumbers = [1:nModes]';

            % concatenate both
            titleNumbers = [mNumbers f]';

            % format titles
            titles = strread(sprintf('Mode %i : f = %8.2f Hz\n', ...
                                     titleNumbers), '%s',        ...
                             'delimiter', '\n');

            % call the plotDisplaced-method of the plot_nodes-class
            self.plot_nodes.plotDisplaced(self.nodes, ...
                                          nodeDisplacements, ...
                                          self.nAdj, titles);

        end

        % ADD NODAL BOUNDARY CONDITION
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        %   p_bc     - bcs for each DOF (1 = clamped, 0 = free)
        %              [bc_u bc_w bc_beta]
        function self = addNodeBC(self, p_i_node, p_bc)

            self.bc(p_i_node, :) = p_bc;

        end

        % ADD NODAL BOUNDARY CONDITION FOR INFINITE ELEMENTS
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        %   p_bc     - spring stiffness and damping coefficient
        %              [Kinf Dinf]
        function self = addNodeBC_inf(self, p_i_node, p_bc)

            self.bc_inf(p_i_node, :) = p_bc;

        end

        % CLAMPED NODAL BC (ALL 3 DOFS FIXED)
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        function self = nodeBC_clamped(self, p_i_node)

            self.addNodeBC(p_i_node, [1 1 1]);

        end

        % JOINTED NODAL BC (ROTATIONAL DOF FREE)
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        function self = nodeBC_jointed(self, p_i_node)

            self.addNodeBC(p_i_node, [1 1 0]);

        end

        % RETURN VECTOR OF FIXED DOFS
        function fixedDOFs = getFixedDOFs(self)

            % 1 -> DOF is fixed
            % 0 -> DOF is free

            % assemble fixed dof vector by taking the transpose of the
            % bc matrix ...
            fixedDOFs = self.bc';
            % ... and concatenate all columns
            fixedDOFs = fixedDOFs(:);

        end

        % ADD NODAL FORCES AND MOMENTS
        %
        % Inputs:
        %   p_i_node - node index the force is acting on
        %   p_f      - forces for each DOF
        %              [f_x f_z m_xz]
        function self = addNodeForces(self, p_i_node, p_f)

            self.fh(p_i_node,:) = p_f;

        end

        % BUILD SYSTEM FORCE VECTOR
        function self = buildFSys(self)

            % build system force vector by reshaping the harmonic force
            % matrix fh in a way that it becomes a column vector with
            % sysDOF elements
            self.fSys = reshape(self.fh.',1,[]).';

        end

    end

end
