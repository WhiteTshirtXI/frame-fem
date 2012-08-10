classdef c_fem_el_beams < handle
%C_FEM_EL_BEAMS - Class containing beam elements of the fem system
% This class contains properties and functions for the handling of the
% FEM system beam-elments. The beam elements are shear-stiff and carry
% normal, transverse and bending moment loads and may have piezoelectric
% properties in the direction of the beam axis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_fem_el_beams.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Properties :
%    n1   - vector of first element node indices
%    n2   - vector of second element node indices
%
%    l    - vector of element lengths
%    a    - vector of element orientation angles
%    rhoA - vector of element mass loading
%    EA   - vector of element axial stiffnesses
%    EI   - vector of element bending stiffnesses
%    dp   - vector of element piezo-electric couling constants
%    epsA - vector of element area-permittivity constants
%
%    ile  - last added element index
%
% Methods :
%    c_fem_el_beams    - constructor
%
%    addBeamPiezo      - add beam(s)
%    addBeamPiezoConst - add beam(s) with constant properties between
%                        nodes
%    addBeam           - add beam(s) w/o piezoelectric coupling
%    addBeamConst      - add beams with constant properties between
%                        nodes w/o piezoelectric coupling
%
%
%    mEl               - return element mass matrix
%    kEl               - return element stiffness matrix
%    cdEl              - return element displacement matrix for beam
%                        element centroid displacements
%    sdEl              - return element stress-diplacement matrix for
%                        beam element centroid stresses
%    mSys              - return system matrices
%
%    elCtrDisp         - return element centroid displacements
%    elCtrStress       - recover element centroid stresses
%    power             - calculate vibrational power in elements
%
%    adj               - return adjacency matrix of inter-connected nodes
%    maxIdx            - return highest node index
%
%
% Static Methods :
%    transGlo       - transform element matrix to global coordinates
%    transGloNode   - transform element node matrix to global
%                     coordinates
%    sysIdx_el      - return system dof indices for node(s) in element
%                     node formulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-06-22 13:10 CEST
% Last Modified : 2012-08-10 10:05 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)
        
        % vector of first element node indices
        n1;
        % vector of second element node indices
        n2;

        % vector of element lengths
        l;
        % vector of element orientation angles
        a;
        % vector of element mass loading
        rhoA;
        % vector of element axial stiffnesses
        EA;
        % vector of element bending stiffnesses
        EI;
        % vector of element piezo-electric coupling constants
        dp;
        % vector of element area-permittivity constants
        epsA;

        % last added element index
        ile;

    end

    % METHODS
    methods
        
        % CONSTRUCTOR
        function s = c_fem_el_beams( p_nElBeams )

            s.n1 = sparse(p_nElBeams, 1);
            s.n2 = sparse(p_nElBeams, 1);
            s.l = sparse(p_nElBeams, 1);
            s.a = sparse(p_nElBeams, 1);
            s.rhoA = sparse(p_nElBeams, 1);
            s.EA = sparse(p_nElBeams, 1);
            s.EI = sparse(p_nElBeams, 1);
            s.dp = sparse(p_nElBeams, 1);
            s.epsA = sparse(p_nElBeams, 1);

            % last added element index = 0
            s.ile = 0;

        end

        % ADD BEAM(S)
        %
        % Inputs:
        %   p_n1   - first node index
        %   p_n2   - second node index
        %   p_l    - element length
        %   p_a    - element orientation angle 
        %   p_rhoA - mass loading
        %   p_EA   - axial stiffness
        %   p_EI   - bending stiffness
        %   p_dp   - piezoelectric coupling factor
        %   p_epsA - area-permittivity constant
        function s = addBeamPiezo( s, p_n1, p_n2, p_l, p_a, p_rhoA,  ...
                                              p_EA, p_EI, p_dp, p_epsA )

            % number of new elements
            nNew = numel(p_n1);

            % new elements index range
            eli = s.ile+1:s.ile+nNew;

            % add element properties to the corresponding vectors
            s.n1(eli) = p_n1(:);
            s.n2(eli) = p_n2(:);
            s.l(eli) = p_l(:);
            s.a(eli) = p_a(:);
            s.rhoA(eli) = p_rhoA(:);
            s.EA(eli) = p_EA(:);
            s.EI(eli) = p_EI(:);
            s.dp(eli) = p_dp(:);
            s.epsA(eli) = p_epsA(:);
            
            % update index of last added element
            s.ile = s.ile+nNew;

        end

        % ADD BEAMS WITH CONSTANT PROPERTIES BETWEEN NODES
        %
        % Inputs:
        %   p_n    - node indices
        %   p_l    - element length (scalar)
        %   p_a    - element orientation angle (scalar)
        %   p_rhoA - mass loading (scalar)
        %   p_EA   - axial stiffness (scalar)
        %   p_EI   - bending stiffness (scalar)
        %   p_dp   - piezoelectric coupling factor (scalar)
        %   p_epsA - area-permittivity constant (scalar)
        function s = addBeamPiezoConst( s, p_n, p_l, p_a, p_rhoA,    ...
                                              p_EA, p_EI, p_dp, p_epsA )

            % one-vector for the constant properites
            one = ones(1,numel(p_n)-1);

            % add elements using artificial properties vectors
            s.addBeamPiezo(p_n(1:end-1), ... % first node index vector
                           p_n(2:end),   ... % second node index vector
                           one*p_l,      ...
                           one*p_a,      ...
                           one*p_rhoA,   ...
                           one*p_EA,     ...
                           one*p_EI,     ...
                           one*p_dp,     ...
                           one*p_epsA);

        end


        % ADD BEAM(S) W/O PIEZOELECTRIC COUPLING
        %
        % Inputs:
        %   p_n1   - first node index
        %   p_n2   - second node index
        %   p_l    - element length
        %   p_a    - element orientation angle 
        %   p_rhoA - mass loading
        %   p_EA   - axial stiffness
        %   p_EI   - bending stiffness
        function s = addBeam( s, p_n1, p_n2, p_l, p_a, p_rhoA, p_EA, ...
                                                                 p_EI )

            % zero-vector for non-effective piezo-electric properites
            zero = zeros(1,numel(p_n1));

            % add elements with zero piezo-electric prpoerties
            s.addBeamPiezo( p_n1, p_n2, p_l, p_a, p_rhoA, p_EA,     ...
                                                    p_EI, zero, zero );

        end

        % ADD BEAMS WITH CONSTANT PROPERTIES BETWEEN NODES W/O
        % PIEZOELECTRIC COUPLING
        %
        % Inputs:
        %   p_n    - node indices
        %   p_l    - element length (scalar)
        %   p_a    - element orientation angle (scalar)
        %   p_rhoA - mass loading (scalar)
        %   p_EA   - axial stiffness (scalar)
        %   p_EI   - bending stiffness (scalar)
        function s = addBeamConst( s, p_n, p_l, p_a, p_rhoA, p_EA,   ...
                                                                  p_EI )

            s.addBeamPiezoConst(p_n, p_l, p_a, p_rhoA, p_EA, p_EI,   ...
                                                                  0, 0);

        end

        % RETURN ELEMENT MASS MATRIX
        %
        % Inputs:
        %   p_el - element index
        function m = mEl( s, p_el )

            m = sparse(8,8);

            rhoA = s.rhoA(p_el);
            l = s.l(p_el);
            ll = l^2;
            rhoAlr420 = rhoA*l/420;

            m(1,1) = rhoAlr420*140;
            m(1,2) = rhoAlr420* 70;
            
            m(2,2) = rhoAlr420*140;

            m(3,3) = rhoAlr420*156;
            m(3,4) = rhoAlr420* 54;
            m(3,5) = rhoAlr420*(-22*l);
            m(3,6) = rhoAlr420*13*l;

            m(4,4) = rhoAlr420*156;
            m(4,5) = rhoAlr420*(-13*l);
            m(4,6) = rhoAlr420*22*l;
            
            m(5,5) = rhoAlr420*4*ll;
            m(5,6) = rhoAlr420*(-3*ll);

            m(6,6) = rhoAlr420*4*ll;

            % make matrix symmetric
            m(2,1) = m(1,2);
            
            m(4,3) = m(3,4);
            m(5,3) = m(3,5);
            m(6,3) = m(3,6);

            m(5,4) = m(4,5);
            m(6,4) = m(4,6);

            m(6,5) = m(5,6);

        end

        % RETURN ELEMENT STIFFNESS MATRIX
        %
        % Inputs:
        %   p_el - element index
        function k = kEl( s, p_el )

            k = sparse(8,8);

            EA = s.EA(p_el);
            EI = s.EI(p_el);
            rl = 1/s.l(p_el);
            EArl = EA*rl;
            EIrl = EI*rl;                     
            EIrll = EIrl*rl;
            EIrlll = EIrll*rl;
            EAdprl = EArl*s.dp(p_el);
            epsArl = rl*s.epsA(p_el);

            k(1,1) =  EArl;
            k(1,2) = -EArl;
            k(1,7) =  EAdprl;
            k(1,8) = -EAdprl;

            k(2,2) =  EArl;
            k(2,7) = -EAdprl;
            k(2,8) =  EAdprl;

            k(3,3) =  12*EIrlll;
            k(3,4) = -12*EIrlll;
            k(3,5) = - 6*EIrll;
            k(3,6) = - 6*EIrll;

            k(4,4) = 12*EIrlll;
            k(4,5) =  6*EIrll;
            k(4,6) =  6*EIrll;
            
            k(5,5) = 4*EIrl;
            k(5,6) = 2*EIrl;

            k(6,6) = 4*EIrl;

            k(7,7) = -epsArl;
            k(7,8) =  epsArl;

            k(8,8) = -epsArl;

            % make matrix symmetric
            k(2,1) = k(1,2);
            k(7,1) = k(1,7);
            k(8,1) = k(1,8);

            k(7,2) = k(2,7);
            k(8,2) = k(2,8);

            k(4,3) = k(3,4);
            k(5,3) = k(3,5);
            k(6,3) = k(3,6);

            k(5,4) = k(4,5);
            k(6,4) = k(4,6);
            
            k(6,5) = k(5,6);

            k(8,7) = k(7,8);
        
        end

        % RETURN ELEMENT DISPLACEMENT MATRIX FOR BEAM ELEMENT CENTROID
        % DISPLACEMENTS
        %
        % Inputs:
        %   p_el - element index
        function cdm = cdEl( s, p_el )

            l = s.l(p_el);
            rl = 1/l;
                     
            cdm = [ 4 4     0      0  0  0 0 0 ;
                    0 0     4      4 -l  l 0 0 ;
                    0 0 12*rl -12*rl -2 -2 0 0 ] / 8;
        
        end

        % RETURN ELEMENT STRESS-DISPLACEMENT MATRIX FOR BEAM ELEMENT
        % CENTROID STRESSES
        %
        % Inputs:
        %   p_el - element index
        function sd = sdEl( s, p_el )

            EA = s.EA(p_el);
            EI = s.EI(p_el);
            rl = 1/s.l(p_el);
            EIrl = EI*rl;                     
            EIrll = EIrl*rl;
                     
            % TODO: check for piezo?
            sd = rl*[ -EA EA         0        0      0      0 0 0 ;
                        0  0 -12*EIrll 12*EIrll 6*EIrl 6*EIrl 0 0 ;
                        0  0         0        0    -EI     EI 0 0 ];
        
        end

        % RETURN SYSTEM MATRICES
        function [M,K] = mSys( s )

            % calculate maximum system matrix size using the highest
            % node index
            maxIdx = s.maxIdx()*4;

            % pre-allocate result matrices
            M = sparse(maxIdx,maxIdx);
            K = sparse(maxIdx,maxIdx);

            % run through all elements
            for eli = 1:s.ile

                % element node indices
                idx = s.sysIdx_el([ s.n1(eli) s.n2(eli) ]);

                % add element matrices to system matrices
                M(idx,idx) = M(idx,idx) + ...
                             s.transGlo(s.mEl(eli), s.a(eli));
                K(idx,idx) = K(idx,idx) + ...
                             s.transGlo(s.kEl(eli), s.a(eli));

            end % elements loop

        end

        % RETURN ELEMENT CENTROID DISPLACEMENTS
        %
        % Inputs:
        %   p_e - element indices
        %   p_u - nodal displacements
        function uc = elCtrDisp( s, p_e, p_u )

            % pre-allocate displacement vector
            uc = zeros(3,numel(p_e));

            % initialize local displacement vector
            u = zeros(8,1);

            % loop through element indices
            for el = 1:numel(p_e)

                % element index
                eli = p_e(el);

                % assign node displacements
                u(1:2:end) = p_u(:,s.n1(eli));
                u(2:2:end) = p_u(:,s.n2(eli));

                % calculate centroid displacement using
                % centroid-displacement-matrix
                uc(:,el) = s.transGloNode(s.cdEl(eli), s.a(eli))*u;

            end % element index loop

        end

        % RECOVER ELEMENT CENTROID STRESSES
        %
        % Inputs:
        %   p_e - element indices
        %   p_u - nodal displacements
        function sig = elCtrStress( s, p_e, p_u )

            % pre-allocate stress vector
            sig = zeros(3,numel(p_e));

            % initialize local displacement vector
            u = zeros(8,1);

            % loop through element indices
            for el = 1:numel(p_e)

                % element index
                eli = p_e(el);

                % assign node displacements
                u(1:2:end) = p_u(:,s.n1(eli));
                u(2:2:end) = p_u(:,s.n2(eli));

                % calculate centroid stresses using
                % stress-displacement-matrix
                sig(:,el) = s.transGloNode(s.sdEl(eli), s.a(eli))*u;

            end % element index loop

        end

        % CALCULATE VIBRATIONAL POWER IN ELEMENTS
        %
        % Inputs:
        %   p_e - element indices
        %   p_u - nodal displacement amplitudes
        %   p_o - excitation frequency
        function P = power( s, p_e, p_u, p_o )

            % calculate stress and velocity amplitudes
            sig = s.elCtrStress( p_e, p_u );
            v = 1i*p_o*s.elCtrDisp( p_e, p_u );

            % calculate powers
            P = 0.5*sum(real(sig.*conj(v)),1);

        end

        % RETURN ADJACENCY MATRIX OF INTER-CONNECTED NODES
        %
        % Inputs:
        %   p_n - node indices, OPTIONAL
        function A = adj( s, p_n )

            % calculate maximum adjacency matrix size using the
            % highest node index
            maxIdx = s.maxIdx();

            if nargin < 2

                p_n = 1:maxIdx;

            end
                
            % pre-allocate adjacency matrix
            A = sparse(maxIdx,maxIdx);

            % fetch adjacency indices using the start and end nodes of
            % each beam
            aIdx = sub2ind([maxIdx,maxIdx], s.n1(1:s.ile),           ...
                                            s.n2(1:s.ile));

            A(aIdx) = 1:s.ile;

            % copy upper triangle to lower triangle
            A = A + A.';

            % return only the rows equivalent to nodes with indices p_n
            A = A(p_n,:);

        end

        % RETURN HIGHEST NODE INDEX
        function idx = maxIdx( s )

            idx = max([ s.n1(1:s.ile) ; s.n2(1:s.ile) ]);

        end

    end

    % STATIC METHODS
    methods (Static)
        

        % TRANSFORM ELEMENT MATRIX TO GLOBAL COORDINATES
        %
        % Inputs:
        %   p_M - element matrix
        %   p_a - element orientation angle
        function M = transGlo( p_M, p_a )

            sa = sin(p_a);
            ca = cos(p_a);

            % transformation matrix
            Te = [ ca  0 -sa   0 0 0 0 0 ;
                    0 ca   0 -sa 0 0 0 0 ;
                   sa  0  ca   0 0 0 0 0 ;
                    0 sa   0  ca 0 0 0 0 ;
                    0  0   0   0 1 0 0 0 ;
                    0  0   0   0 0 1 0 0 ;
                    0  0   0   0 0 0 1 0 ;
                    0  0   0   0 0 0 0 1 ];

            % transpose of transformation matrix
            TeT = [  ca   0 sa  0 0 0 0 0 ;
                      0  ca  0 sa 0 0 0 0 ;
                    -sa   0 ca  0 0 0 0 0 ;
                      0 -sa  0 ca 0 0 0 0 ;
                      0   0  0  0 1 0 0 0 ;
                      0   0  0  0 0 1 0 0 ;
                      0   0  0  0 0 0 1 0 ;
                      0   0  0  0 0 0 0 1 ];

            % transform matrix
            M = TeT*p_M*Te;

        end

        % TRANSFORM ELEMENT NODE MATRIX TO GLOBAL COORDINATES
        %
        % Inputs:
        %   p_M - element node matrix
        %   p_a - element orientation angle
        function M = transGloNode( p_M, p_a )

            sa = sin(p_a);
            ca = cos(p_a);

            % transformation matrix
            Te = [ ca  0 -sa   0 0 0 0 0 ;
                    0 ca   0 -sa 0 0 0 0 ;
                   sa  0  ca   0 0 0 0 0 ;
                    0 sa   0  ca 0 0 0 0 ;
                    0  0   0   0 1 0 0 0 ;
                    0  0   0   0 0 1 0 0 ;
                    0  0   0   0 0 0 1 0 ;
                    0  0   0   0 0 0 0 1 ];

            % transpose of transformation matrix
            TeT = [  ca sa 0 ;
                    -sa ca 0 ;
                     0   0 1 ];

            % transform matrix
            M = TeT*p_M*Te;

        end

        % RETURN SYSTEM DOF INDICES FOR NODE(S) IN ELEMENT NODE
        % FORMULATION
        %
        % Inputs:
        %   p_i - node index or index vector
        function idx = sysIdx_el( p_i )

            % number of nodal degrees of freedom
            NDOF = 4;
            % local dof index vector
            ldiv = 0:(NDOF-1);
            
            % first: multiply the element node index vector by NDOF and
            %        substract (NDOF-1) to get the index vector in terms
            %        of global DOF indices
            idx_g = NDOF*p_i(:).'-(NDOF-1);
            
            % second: expand the node DOF index row vector idx_g to get
            %         a matrix with NDOF rows
            idx = repmat(idx_g, NDOF, 1);

            % third: expand the local DOF index vector [0;1;2] in the
            %        same way and sum up both matrices
            idx = idx + repmat(ldiv', 1, numel(p_i));

            % fourth: the final dof index vector is generated by
            %         concatenating all of the matrix' columns
            idx = idx(:);

            % re-shape index vector to index matrix
            idx = reshape(idx, NDOF, []).';

            % concatenate all columns
            idx = idx(:);

        end

    end


end
