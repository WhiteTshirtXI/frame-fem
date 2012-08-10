classdef c_vv_beam_sInf < handle
%C_VV_BEAM_SINF - Verification class: semi-infinite beam
% This class is for verifying the fem system using the forced vibrations
% of a semi-infinite beam which has exact analytical solutions as a
% reference. The example beam properties are taken from:
% WANG, C. and J.C.S. LAI (2000): Modelling the Vibration Behaviour of
% Infinite Structures By FEM. Journal of Sound and Vibration 229(3),
% pp. 453-466.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_vv_beam_sInf.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Constant Properties :
%    H   - beam cross section height (m)
%    W   - beam cross section width (m)
%    L   - beam length (m)
%    A   - beam angle (rad)
%    RHO - density (kg/m^3)
%    E   - young's modulus (N/m^2)
%    OM  - excitation angular frequency (rad/s)
%    F   - excitation force amplitudes (N)
%    NEL - element density on beam (1/m)
%
% Properties :
%    omega    - excitation frequency
%    l        - length
%
%    rhoA     - mass
%    EA       - axial stiffness
%    EI       - bending stiffness
%
%    nodes    - beam nodes
%    beam     - beam
%    frame    - frame class
%    fem      - fem class
%
%    nInfinit - infinite element nodes
%
% Methods :
%    c_vv_beam_sInf - constructor
%
%    studyRefine    - refinement study
%    studyLength    - beam length study
%    vv             - perform vv-study
%
%    initFrame      - initialize frame class
%    initFEA        - initialize fe-analysis 
%    initNodes      - initialize nodes
%    initBeam       - initialize beam specification
%    initBCs        - initialize boundary conditions
%    initForces     - initialize nodal forces
%
%    dpm            - calculate driving point mobility
%    dpmA           - analytical solution for driving point mobility
%                     matrix
%
%    txtOut         - text output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-07-19 10:46 CEST
% Last Modified : 2012-08-10 14:28 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CONSTANT PROPERTIES %
    properties (Constant)

        % DEFAULT BEAM PROPERTIES

        % square cross section: height and width (m)
        H =   1e-3;
        W =  10e-3;

        % beam length (m)
        L = 200e-3;

        % beam angle (rad)
        A = 0;

        % density (kg/m^3)
        RHO = 7860;

        % young's modulus (N/m^2)
        E = 210000e6;

        % excitation angular frequency (rad/s)
        OM = 100*2*pi;

        % excitation force amplitudes (N)
        F = [0 ; 1 ; 0];

        % DEFAULT FEM PROPERTIES
        
        % element density on beam (1/m)
        NEL = 10;

    end

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)

        % excitation frequency
        omega;

        % length
        l;

        % DERIVED PROPERTIES

        % mass
        rhoA;
        
        % stiffnesses
        EA;
        EI;

        % FEM PROPERTIES

        % beam nodes
        nodes;
        % beam
        beam;
        % frame class
        frame;
        % fem class
        fem;

        % infinite element nodes
        nInfinit;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        function s = c_vv_beam_sInf()

            s.rhoA = s.RHO*s.H*s.W;
            s.EA = s.E*s.H*s.W;
            s.EI = s.E*s.H^3*s.W/12;

            s.l = s.L;
            s.omega = s.OM;

            s.initNodes();
            s.initBeam();

            s.nInfinit = [2];

        end

        % REFINEMENT STUDY
        %
        % Inputs:
        %   p_r - number of refinement levels (OPTIONAL)
        function [dpmRef,dpmD] = studyRefine( s, p_r )

            if ~exist('p_r','var'), p_r = 1 ; end

            % calculate real part of analytical driving point moblity
            % matrix
            dpmRef = real(s.dpmA(s.omega));

            % pre-allocate results matrix (relative difference between
            % numeric and analytical solution)
            dpmD = zeros(3,3,p_r);

            % perfom refinement study
            for r = 1:p_r

                % initialize refined beam
                s.initBeam(r);

                % calculate real part of driving point mobility matrix
                % and compute difference to analytical solution
                dpmD(:,:,r) = max(real(s.dpm()),eps)./max(dpmRef,eps)-1;

            end

            % reset beam and FEA-system
            s.initBeam();
            s.initFEA();

        end

        % BEAM LENGTH STUDY
        %
        % Inputs:
        %   p_l - number of length increments (OPTIONAL)
        %   p_r - number of refinement levels (OPTIONAL)
        function [dpmRef,dpmD] = studyLength( s, p_l, p_r )

            if ~exist('p_l','var'), p_l = 1 ; end
            if ~exist('p_r','var'), p_r = 1 ; end

            % calculate real part of analytical driving point moblity
            % matrix
            dpmRef = real(s.dpmA(s.omega));

            % pre-allocate results matrix (relative difference between
            % numeric and analytical solution)
            dpmD = zeros(3,3,p_r,p_l);

            % perfom length incrementing
            for l = 1:p_l

                % increase beam length
                s.l = s.L*2.^(l-1);

                % initialize nodes
                s.initNodes();

                % perform refinement study
                [~,dpmD(:,:,:,l)] = s.studyRefine(p_r);

            end

            % reset beam length
            s.l = s.L;
            s.initNodes();
            s.initBeam();
            s.initFEA();

        end

        % PERFORM VV-STUDY
        function s = vv( s )

            % VV OPTIONS
            % number of length increments
            nl = 4;

            % number of refinement levels
            nr = 3;

            % calculate number of subplots
            nPltRow = ceil(sqrt(nl));
            nPltCol = ceil(nl/nPltRow);

            % do beam length study
            [ref,d] = s.studyLength( nl, nr );

            % create parent figure
            figure;
            title 'Refinement study - semi infinite beam'

            % crate subfigure for each length increment
            for jl = 1:nl

                % reshape results difference matrix for nicer plotting
                dPlot = abs(reshape(d(:,:,:,jl),9,nr)).';

                % plot grid convergence subfigure
                subplot(nPltRow,nPltCol,jl);
                loglog(2.^[0:nr-1]*s.NEL,dPlot(:,[1 5 6 8 9]),'-*');
                title(sprintf('Beam length %.2f m', s.L*2.^(jl-1)));
                xlabel 'Element density (1/m)';
                ylabel 'Relative error';

            end
            
            legend({'\epsilon(Y_{uN})',...
                   '\epsilon(Y_{wQ})',...
                   '\epsilon(Y_{{\beta}Q})',...
                   '\epsilon(Y_{wM})',...
                   '\epsilon(Y_{{\beta}M})'},...
                   'Position',[0.85 0.4 0.09 0.26])

        end

        % INIALIZE FRAME CLASS
        function s = initFrame(s)

            % initialize frame-definition class
            s.frame = c_frame_def(s.nodes);

            % add the beams to frame-class
            s.frame.addBeam(s.beam);

        end

        % INITIALIZE FE-ANALYSIS
        function s = initFEA(s)

            s.initFrame();
            s.initBCs();
            s.initForces();
            s.fem = s.frame.discretize();

        end

        % INITIALIZE NODES
        %
        % Inputs:
        %   p_a - beam angle (OPTIONAL)
        function s = initNodes( s, p_a )

            if nargin < 2
                p_a = s.A;
            end

            s.nodes = [0 0 ; cos(p_a) sin(p_a)]*s.l;

        end

        % INITIALIZE BEAM SPECIFICATION
        %
        % Inputs:
        %   p_r - refinement level (OPTIONAL)
        function s = initBeam( s, p_r )

            if nargin < 2
                p_r = 1;
            end
            
            s.beam = [ 1 2 s.rhoA s.EA s.EI s.NEL*2^(p_r-1) ];

        end

        % INITIALIZE BOUNDARY CONDITIONS
        %
        % Inputs:
        %   p_o - excitation frequency (OPTIONAL)
        %
        function s = initBCs(s, p_o)

            if nargin < 2
                p_o = s.OM;
            end

            s.omega = p_o;

            % initialize infinite element boundaries
            s.frame.nodeBC_infinite(s.nInfinit, s.omega);

        end

        % INITIALIZE NODAL FORCES
        %
        % Inputs:
        %   p_f - force vector (OPTIONAL)
        function s = initForces(s, p_f)

            if nargin < 2
                p_f = s.F;
            end

            s.frame.addHarmonicForce(1, p_f);

        end

        % CALCULATE DRIVING POINT MOBILITY
        function dpm = dpm( s )

            % initialize driving point mobility matrix
            dpm = zeros(3);

            % run through three different load cases where a single unit
            % force is applied at the first beam node cycling through
            % normal, transversal and moment loading.
            
            jf = 0;

            for f = eye(3)

                jf = jf + 1;

                s.initFrame();
                s.initBCs();
                s.initForces(f);
                s.fem = s.frame.discretize();

                % calculate complex velocity amplitudes
                v = 1i*s.omega*s.fem.harmonicAnalysis(s.omega);

                % calculate driving point mobilities (driving point =
                % first node)
                dpm(:,jf) = v(1:3,1);

            end

        end


        % ANALYTICAL SOLUTION FOR DRIVING POINT MOBILITY MATRIX
        %
        % Inputs:
        %   p_o - excitation frequency
        function dpm = dpmA( s, p_o )

            % wave numbers
            kL = p_o*sqrt(s.rhoA/s.EA);
            kB = sqrt(p_o)*(s.rhoA/s.EI)^.25;

            % driving point mobility matrix according to beam theory
            dpm = [ 1./(s.EA*kL) 0 0  ;
                    0 (1-1i)./(s.EI*kB.^3) (1-2i)./(s.EI*kB.^2) ;
                    0 1./(s.EI*kB.^2)      (1-1i)./(s.EI*kB) ].*p_o;

        end

        % TEXT OUTPUT
        %
        % Inputs:
        %   p_t - text string
        function s = txtOut(s, p_t)

            fprintf(p_t);

        end
        
    end

end
