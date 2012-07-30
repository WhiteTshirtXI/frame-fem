classdef c_vv_beam_nModes < handle
%C_VV_BEAM_NMODES - Verification class: beam normal modes
% This class is for verifying the fem system using the free vibrations
% of a single beam which has exact analytical solutions as a reference.
% The example beam properties are taken from:
% FAHY (2007): Sound and Structural Vibration: Radiation, Transmission
% and Response. Elsevier. Pp. 460-469.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_vv_beam_nModes.m
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
%    NEL - element density on beam (1/m)
%    NM  - default maximum modes number
%
% Properties :
%    rhoA     - mass
%    EA       - axial stiffness
%    EI       - bending stiffness
%
%    nodes    - beam nodes
%    beam     - beam
%    frame    - frame class
%    fem      - fem class
%
%    nSpprted - simply supported nodes
%    nClamped - clamped nodes
%
% Methods :
%    c_vv_beam_nModes - constructor
%
%    studyRefine      - refinement study
%    studyAlpha       - beam angle study
%    vv               - perform vv-study
%
%    initFrame        - initialize frame class
%    initFEA          - initialize fe-analysis 
%    initNodes        - initialize nodes
%    initBeam         - initialize beam specification
%    initBCs          - initialize boundary conditions
%
%    eigFa            - analytical solution for eigenfrequencies
%
%    txtOut           - text output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-07-03 15:00 CEST
% Last Modified : 2012-07-30 16:40 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CONSTANT PROPERTIES %
    properties (Constant)

        % DEFAULT BEAM PROPERTIES

        % square cross section: height and width (m)
        H =   5e-3;
        W =  20e-3;

        % beam length (m)
        L = 500e-3;

        % beam angle (rad)
        A = 0;

        % density (kg/m^3)
        RHO = 2700;

        % young's modulus (N/m^2)
        E = 71000e6;

        % DEFAULT FEM PROPERTIES
        
        % element density on beam (1/m)
        NEL = 10;

        % default maximum modes number
        NM = 4;

    end

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)

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

        % simply supported nodes
        nSpprted;
        % clamped nodes
        nClamped;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        function s = c_vv_beam_nModes()

            s.rhoA = s.RHO*s.H*s.W;
            s.EA = s.E*s.H*s.W;
            s.EI = s.E*s.H^3*s.W/12;

            s.initNodes();
            s.initBeam();

            s.nSpprted = [2];
            s.nClamped = [1];

        end

        % REFINEMENT STUDY
        %
        % Inputs:
        %   p_n - maximum number of modes (OPTIONAL)
        %   p_r - number of refinement levels (OPTIONAL)
        function [f,d] = studyRefine( s, p_n, p_r )

            if ~exist('p_n','var'), p_n = s.NM ; end
            if ~exist('p_r','var'), p_r = 1 ; end

            % maximum node number
            nMax = p_n;

            % calculate analytical eigenfrequency vector
            f = s.eigFa(nMax);

            % pre-allocate results matrix (relative difference between
            % numeric and analytical solution)
            d = zeros(p_n,p_r);

            % perfom refinement study
            for r = 1:p_r

                % initialize refined beam and FEA-system
                s.initBeam(r);
                s.initFEA();

                % if neccessary, correct maximum node number
                nMax = min(nMax, s.fem.maxModes());

                % calculate eigenfrequencies of discretized system and
                % compute difference to analytical solution
                d(1:nMax,r) = max(s.fem.eigF(nMax),eps)./ ...
                              max(f(1:nMax),eps)-1;

            end

            % reset beam and FEA-system
            s.initBeam();
            s.initFEA();

            % return results
            f = f(1:nMax);
            d = d(1:nMax,:);

        end

        % BEAM ANGLE STUDY
        %
        % Inputs:
        %   p_a - angle list
        %   p_n - maximum number of modes (OPTIONAL)
        %   p_r - number of refinement levels (OPTIONAL)
        function [f,d] = studyAlpha( s, p_a, p_n, p_r )

            if ~exist('p_a','var'), p_a = s.A ; end
            if ~exist('p_n','var'), p_n = s.NM ; end
            if ~exist('p_r','var'), p_r = 1 ; end

            % pre-allocate results matrix (relative difference between
            % numeric and analytical solution)
            d = zeros(p_n,p_r,numel(p_a));

            % run through all specifiec angles
            for ia = 1:numel(p_a)

                a = p_a(ia);
            
                % change nodes
                s.initNodes(a);

                % perform refinement study
                [f,dTmp] = s.studyRefine(p_n, p_r);

                d(1:numel(f),:,ia) = dTmp;

                ia = ia + 1;

            end

            d = d(1:numel(f),:,:);

        end

        % PERFORM VV-STUDY
        function s = vv( s )

            % VV OPTIONS
            % select random angles
            alpha = rand(10,1)*2*pi;

            % number of modes
            nm = 8;

            % number of refinement levels
            nr = 6;

            [f,d] = s.studyAlpha( alpha, nm, nr );

            % check if the standard deviation in the results exceeds a
            % specific tolerance
            dev = max(max(std(d,0,3)));
            devTol = sqrt(eps);

            if dev > devTol

                s.txtOut(['*** WARNING : results not independent '...
                          'of beam angle!\n  -> dev = '...
                          sprintf('%8.2e > %8.2e', [dev devTol])]);

            else

                s.txtOut(['Max. standard deviation in beam angle '...
                          'study: '...
                          sprintf('%8.2e < %8.2e', [dev devTol]) ...
                          '\n  -> OK!\n']);

                % plot grid convergence figure
                figure;
                loglog(2.^[0:nr-1]*s.NEL,abs(mean(d,3)).','-*');
                title 'Refinement study - eigenfrequencies of a beam'
                xlabel 'Number of elements'
                ylabel 'Error'
                legend(num2str(f,'%7.1f Hz'));

            end

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

            s.nodes = [0 0 ; cos(p_a) sin(p_a)]*s.L;

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
        function s = initBCs(s)

            % initialize simply supported boundaries
            s.frame.nodeBC_spprted(s.nSpprted);

            % initialize clamped boundaries
            s.frame.nodeBC_clamped(s.nClamped);

        end


        % ANALYTICAL SOLUTION FOR EIGENFREQUENCIES
        %
        % Inputs:
        %   p_n - number of desired eigenfrequencies
        function f = eigFa( s, p_n )

            % REFERENCE:
            % YOUNG, D. and R. FELGAR: Tables of characteristic
            % functions representing normal modes of vibration of a
            % beam. Engineering Research Series. University of Texas,
            % 1949.
            
            % pre-allocate frequency vector
            f = zeros(p_n,1);

            % get longitudinal modes eigenfrequencies according to
            % boundary conditions
            if numel(s.nSpprted) + numel(s.nClamped) == 1
                
                % CASE 1 : one-sidedly fixed rod
                f_l = inline('(2*(1:n).''-1)*pi/2');

            else

                % CASE 2 : both sides fixed or free
                f_l = inline('(1:n).''*pi');

            end

            % get bending modes eigenfrequencies according to boundary
            % conditions
            switch numel(s.nSpprted)

                % CASE 1 : supported-supported beam
                case 2

                    % inline function for n bending eigenfrequencies
                    f_b = inline('((1:n).''*pi).^2');

                % CASE 2 : clamped-supported or free-supported beam
                case 1

                    f_b = inline(['[  3.92660230 ;' ...
                                  '   7.06858275 ;' ...
                                  '  10.21017613 ;' ...
                                  '  13.35176878 ;' ...
                                  '  16.49336143 ;' ...
                                  '  ((4*(6:n).''+1)*pi/4) ].^2']);

                % no simply supported bc
                case 0

                    switch numel(s.nClamped)

                        % CASE 3 : clamped-free beam
                        case 1

                            f_b = inline(['[  1.87510410 ;' ...
                                          '   4.69409113 ;' ...
                                          '   7.85475743 ;' ...
                                          '  10.99554074 ;' ...
                                          '  14.13716839 ;' ...
                                          ' ((2*(6:n).''-1)*pi/2)' ...
                                          ' ].^2']);

                        % CASE 4 : clamped-clamped or free-free beam
                        otherwise

                            f_b = inline(['[  4.73004080 ;' ...
                                          '   7.85320460 ;' ...
                                          '  10.99560730 ;' ...
                                          '  14.13716550 ;' ...
                                          '  17.27875960 ;' ...
                                          '  ((2*(6:n).''+1)*pi/2)' ...
                                          ' ].^2']);

                    end

            end

            % rigid body translations (f = 0)
            f_rb = [];
            if numel(s.nSpprted) == 1 & numel(s.nClamped) == 0
                f_rb = zeros(1,1);
            elseif numel(s.nSpprted) == 0 & numel(s.nClamped) == 0
                f_rb = zeros(3,1);
            end

            % create analytical eigenfrequencies vector
            f = sort([ f_rb ;
                       f_l(p_n)*sqrt(s.EA/s.rhoA) ;                  ...
                       f_b(p_n)*sqrt(s.EI/s.rhoA)/s.L ])/(2*pi*s.L);
            f = f(1:p_n);

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
