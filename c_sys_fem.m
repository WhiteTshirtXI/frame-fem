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
%
%    nodes      - nodes
%    el_beams   - beam elements
%
%    bc_inf     - infinite boundary conditions for each node
%
%    eigVal     - matrix of eigenvalues of the system
%    eigVec     - matrix of eigenvectors of the system
%    eigRecalc  - recalculation of eigenvalues required
%
% Methods :
%    c_sys_fem       - constructor
%    sysDOF          - return overall system degrees of freedom
%    addNodes        - add nodes
%    addElBeam       - add beam element(s)
%    addElBeamComst  - add beam elements with constant properties
%                      between nodes
%
%    harmonicAnalysis - do a harmonic analysis in frequency domain
%    eigCalc         - calculate eigenvalues and eigenvectors
%    eigOmega        - calculate angular eigenfrequencies
%    eigF            - calculate eigenfrequencies
%    maxModes        - maximum number of eigenmodes
%
%    addNodeBC       - add nodal boundary condition
%    addNodeBC_inf   - add nodal boundary condition for infinite
%                      elements
%
%    addNodeForces   - add nodal forces and moments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-18 12:50 CEST
% Last Modified : 2012-07-30 17:11 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)
        % number of DOFs per node
        % (equals 3, because, currently, only beam-truss-like elements
        % are supported)
        NDOF = 3;

        % nodes class
        nodes;

        % beam elements class
        el_beams;

        % infinite boundary conditions for each node
        bc_inf;

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
        %   p_n        - (maximum) number of nodes
        %   p_nElBeams - (maximum) number of beam elements
        function s = c_sys_fem( p_n , p_nElBeams )

            % initialize nodes and element classes
            s.nodes = c_fem_nodes( p_n );
            s.el_beams = c_fem_el_beams( p_nElBeams );

            % initialize infinite bcs to zero
            s.bc_inf = sparse(p_n,2*s.NDOF);

            % overall degrees of freedom
            sysDOF = s.sysDOF();

            % allocate memory for the eigenvalue matrices
            s.eigVal = zeros(sysDOF);
            s.eigVec = zeros(sysDOF);

        end

        % RETURN OVERALL SYSTEM DEGREES OF FREEDOM
        function sysDOF = sysDOF(s)

            sysDOF = s.nodes.n()*s.NDOF;

        end

        % ADD NODES
        %
        % Inputs:
        %   p_x - node x-coordinate(s)
        %   p_z - node z-coordinate(s)
        function idx = addNodes( s, p_x, p_z )

            idx = s.nodes.addNode(p_x, p_z);
            
            % a recalculation of eigenvalues is required
            s.eigRecalc = true;

        end

        % ADD BEAM ELEMENT(S)
        %
        % Inputs:
        %   p_n1   - first node index
        %   p_n2   - second node index
        %   p_l    - element length
        %   p_a    - element orientation angle 
        %   p_rhoA - mass loading
        %   p_EA   - axial stiffness
        %   p_EI   - bending stiffness
        function s = addElBeam(s, p_n1, p_n2, p_l, p_a,        ...
                                                     p_rhoA, p_EA, p_EI)

            s.el_beams.addBeam(p_n1, p_n2, p_l, p_a, p_rhoA, p_EA,...
                                                                  p_EI);
            
            % a recalculation of eigenvalues is required
            s.eigRecalc = true;

        end

        % ADD BEAM ELEMENTS WITH CONSTANT PROPERTIES BETWEEN NODES
        %
        % Inputs:
        %   p_n    - node indices
        %   p_l    - element length (scalar)
        %   p_a    - element orientation angle (scalar)
        %   p_rhoA - mass loading (scalar)
        %   p_EA   - axial stiffness (scalar)
        %   p_EI   - bending stiffness (scalar)
        function s = addElBeamConst(s, p_n, p_l, p_a, p_rhoA,  ...
                                                             p_EA, p_EI)

            s.el_beams.addBeamConst(p_n, p_l, p_a, p_rhoA, p_EA,  ...
                                                                  p_EI);
            
            % a recalculation of eigenvalues is required
            s.eigRecalc = true;

        end

        % DO A HARMONIC ANALYSIS IN FREQUENCY DOMAIN
        %
        % Inputs:
        %   p_om - angular frequency of harmonic excitation forces
        function u = harmonicAnalysis(s, p_om)

            % get index vector of free (non-fixed) nodes
            idx = s.nodes.idxFreeDOFs();

            % pre-allocate displacement vector
            u = zeros(3,s.nodes.iln);

            % assemble system matrices
            [MSys,KSys] = s.el_beams.mSys();
            DSys = sparse(s.sysDOF(),s.sysDOF());

            % add frequency dependent infinite boundary conditions to
            % the diagonals of the system matrices
            KSys = KSys + diag(reshape(s.bc_inf(:,1:2:end).',1,[]));
            DSys = DSys + diag(reshape(s.bc_inf(:,2:2:end).',1,[]));

            % calculate inverse transfer matrix
            % (for a system of linear equations Ax=b this equals the
            % matrix A!)
            HSysI = (-p_om^2*MSys+1j*p_om*DSys+KSys);

            % calculate complex displacement amplitudes via x=(A^-1)b
            u(idx) = HSysI(idx,idx)\s.nodes.f(idx);

        end

        % CALCULATE EIGENVALUES AND EIGENVECTORS
        %
        % Inputs:
        %   p_n - number of eigenmodes to be returned.
        %         OPTIONAL: if non-existent, all eigenmodes will be
        %         returned
        function s = eigCalc(s, p_n)
            
            % check if the parameter p_n has been specified
            if ~exist('p_n','var')
                p_n = s.sysDOF()-1;
            end
                
            % get index vector of free (non-fixed) nodes
            idx = s.nodes.idxFreeDOFs();

            % number of free (non-fixed) DOFs
            bcNFree = numel(idx);

            % initialize eigenvector and eigenvalue matrices as zero
            s.eigVec = zeros(bcNFree,p_n);
            s.eigVal = zeros(bcNFree);

            % assemble system matrices
            [MSys,KSys] = s.el_beams.mSys();

            % calculate eigenvaues and eigenvectors close to zero
            [s.eigVec(idx,:),s.eigVal] = eigs(  KSys(idx,idx), ...
                                                     -MSys(idx,idx), ...
                                                      p_n, 1);

        end

        % CALCULATE ANGULAR EIGENFREQUENCIES
        %
        % Inputs:
        %   p_n - number of eigenmodes to be returned.
        %         OPTIONAL: if non-existent, all eigenmodes will be
        %         returned
        function omega = eigOmega(s, p_n)

            % calculate eigenvectors and eigenfrequencies
            s.eigCalc(p_n);
            
            % return vector of angular eigenfrequencies and re-organize
            % eigenvalue and eigenvector matrices
            [omega,idx] = sort(imag(sqrt(diag(s.eigVal))));
            s.eigVal = s.eigVal(idx,idx);
            s.eigVec = s.eigVec(:,idx);

        end

        % CALCULATE EIGENFREQUENCIES
        %
        % Inputs:
        %   p_n - number of eigenmodes to be returned.
        %         OPTIONAL: if non-existent, all eigenmodes will be
        %         returned
        function f = eigF(s, p_n)

            f = s.eigOmega(p_n)./(2*pi);

        end

        % MAXIMUM NUMBER OF EIGENMODES
        function nMax = maxModes( s )

            nMax = numel(s.nodes.idxFreeDOFs());

        end

        % ADD NODAL BOUNDARY CONDITION
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        %   p_bc     - bcs for each DOF (1 = clamped, 0 = free)
        %              [bc_u bc_w bc_beta]
        function s = addNodeBC(s, p_i_node, p_bc)

            s.nodes.fixNode( p_i_node, p_bc.' );

        end

        % ADD NODAL BOUNDARY CONDITION FOR INFINITE ELEMENTS
        %
        % Inputs:
        %   p_i_node - node index the BC is to be applied to
        %   p_bc     - spring stiffness and damping coefficient
        %              for each DOF
        %              [K_u K_w K_beta D_u D_w D_beta]
        function s = addNodeBC_inf(s, p_i_node, p_bc)

            s.bc_inf(p_i_node, :) = p_bc;

        end

        % ADD NODAL FORCES AND MOMENTS
        %
        % Inputs:
        %   p_i_node - node index the force is acting on
        %   p_f      - forces for each DOF
        %              [f_x f_z m_xz]
        function s = addNodeForces(s, p_i_node, p_f)

            s.nodes.setForce( p_i_node, p_f.' );

        end

    end

end
