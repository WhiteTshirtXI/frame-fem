classdef c_prj_dframe < handle
%C_PRJ_DFRAME - Project class for a simple double-frame structure
% This class contains all geometrical and material specifications for
% the vibrational analysis of a double frame structure using the
% beam-fem-classes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_prj_dframe.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class definition
%
% Properties :
%    D_LF        - default distance between force excitation point and
%                  frame structure (m)
%    D_LINF      - default length of the infinite element boundary beams
%                  (m)
%    D_H         - default height of the frame structure (m)
%    D_W         - default width of the frame structure (bolt length)
%                  (m)
%
%    D_NSKIN     - default element density on skin beams (1/m)
%    D_RHOASKIN  - default mass/length on skin beams (kg/m)
%    D_EASKIN    - default axial rigidity on skin beams (N)
%    D_EISKIN    - default bending rigidity on skin beams (Nm^2)
%
%    D_NFRAME    - default element density on frames (1/m)
%    D_RHOAFRAME - default mass/length on frame beams (kg/m)
%    D_EAFRAME   - default axial rigidity on frame beams (N)
%    D_EIFRAME   - default bending rigidity on frame beams (Nm^2)
%
%    D_NBOLT     - default element density on bolt (1/m)
%    D_RHOABOLT  - default mass/length on bolt beams (kg/m)
%    D_EABOLT    - default axial rigidity on bolt beams (N)
%    D_EIBOLT    - default bending rigidity on bolt beams (Nm^2)
%
%    D_OMEGA     - default excitation force frequency (rad/s)
%
%    D_F         - default excitation force amplitudes (N)
%    D_FNODES    - default excitation force node index
%
%    D_BC_CLAMPED - default clamped (fixed) nodes
%    D_BC_SPPRTED - default simply supported nodes
%    D_BC_INFINIT - default infinite element nodes for harmonic analysis
%
%    D_PROBES     - default probe node indices
%
%    D_CFG_BASHM  - default bashmode configuration setting
%
%    lF          - distance between force excitation point and
%                  frame structure (m)
%    lInf        - length of the infinite element boundary beams
%                  (m)
%    h           - height of the frame structure (m)
%    w           - width of the frame structure (bolt length) (m)
%        
%    nSkin       - element density on skin beams (1/m)
%    rhoASkin    - mass/length on skin beams (kg/m)
%    EASkin      - axial rigidity on skin beams (N)
%    EISkin      - bending rigidity on skin beams (Nm^2)
%        
%    nFrame      - element density on frames (1/m)
%    rhoAFrame   - mass/length on frame beams (kg/m)
%    EAFrame     - axial rigidity on frame beams (N)
%    EIFrame     - bending rigidity on frame beams (Nm^2)
%
%    nBolt       - element density on bolt (1/m)
%    rhoABolt    - mass/length on bolt beams (kg/m)
%    EABolt      - axial rigidity on bolt beams (N)
%    EIBolt      - bending rigidity on bolt beams (Nm^2)
%
%    omega       - excitation force frequency (rad/s)
%
%    f           - excitation force amplitudes (N)
%    fNodes      - excitation force node index
%
%    nodes       - nodes coordinate vector
%    beams       - beams specification vector
%    frame       - frame class
%    fem         - fem-system class
%
%    bc_clamped   - list of clamped (fixed) nodes
%    bc_spprted   - list of simply supported nodes
%    bc_infinit   - list of infinite element nodes for harmonic analysis
%
%    probes       - probe node indices
%
%    cfg_bashm    - bash (quiet) mode activation 
%
% Methods :
%    c_prj_dframe - constructor
%    setDefaults  - set default values
%
%    initFrame    - initialize frame class
%    initFEA      - initialize FE-analysis
%    initNodes    - initialize structural nodes
%    initBeams    - initialize structural beams
%    initBCs      - initialize boundary conditions
%    initForces   - initialize nodal forces
%
%    calcTL       - calculate transmission loss in double-frame
%
%    study_lInf   - parameter study of infinite element beam length
%    study_EIBolt - parameter study of bolt bending rigidity
%
%    getIndices   - return node indices of specific nodes
%    extractNodes - extract specific node variable data from the
%                   discretized system
%    vectIndices  - return data vector indices for specified nodes
%
%    txtOut       - text output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-06-06 16:02 CEST
% Last Modified : 2012-06-12 14:05 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CONSTANT PROPERTIES %
    properties (Constant)

        %%% DEFAULT VALUES %%%

        % lengths
        D_LF = 3.0;
        D_LINF = 5.0;
        D_H = 0.12;
        D_W = 0.1;

        % default skin beam properties
        % element density (elements per unit length)
        D_NSKIN = 10;
        % mass per unit length
        D_RHOASKIN = 5.54;
        % axial rigidity
        D_EASKIN = 143.7e6;
        % bending rigidity
        D_EISKIN = 12266.5;

        % default frame beam properties
        % element density (elements per unit length)
        D_NFRAME = 15;
        % mass per unit length
        D_RHOAFRAME = 1.17;
        % axial rigidity
        D_EAFRAME = 30.24e6;
        % bending rigidity
        D_EIFRAME = 66949.1;

        % default bolt beam properties
        % element density (elements per unit length)
        D_NBOLT = 15;
        % mass per unit length
        D_RHOABOLT = 0.788;
        % axial rigidity
        D_EABOLT = 20.5e6;
        % bending rigidity
        D_EIBOLT = 247.397;

        % frequencies
        D_OMEGA = 100*2*pi;

        % force amplitudes [ Fx ; Fz ; My ]
        D_F = [ 0 ; 100.0 ; 0 ];
        % force node index
        D_FNODES = [ 2 ];

        % boundary conditions
        % clamped (fixed) nodes
        D_BC_CLAMPED = 'none';
        % simply supported nodes
        D_BC_SPPRTED = 'none';
        % infinite element nodes for harmonic analysis
        D_BC_INFINIT = 'firstAndLast';

        % probe node indices (after double-frame structure)
        % TODO: correct to [ 6 ]
        D_PROBES = 6;

        % default configuration options
        % bash mode
        D_CFG_BASHM = 1;

    end

    % PRIVATE PROPERTIES %
    properties (SetAccess = private)

        %%% STRUCTURAL PROPERTIES %%%

        % lengths
        lF;
        lInf;
        h;
        w;

        % skin beam properties
        % element density (elements per unit length)
        nSkin;
        % mass per unit length
        rhoASkin;
        % axial rigidity
        EASkin;
        % bending rigidity
        EISkin;

        % frame beam properties
        % element density (elements per unit length)
        nFrame;
        % mass per unit length
        rhoAFrame;
        % axial rigidity
        EAFrame;
        % bending rigidity
        EIFrame;

        % bolt beam properties
        % element density (elements per unit length)
        nBolt;
        % mass per unit length
        rhoABolt;
        % axial rigidity
        EABolt;
        % bending rigidity
        EIBolt;

        % frequencies
        omega;

        % force amplitudes
        f;
        % force node index
        fNodes;

        %%% STRUCTURAL REPRESENTATION VIA NODES AND BEAMS %%%

        % nodes coordinate vector
        nodes;

        % beams specification vector
        % (format: [n1 n2 rhoA EA EI nEl])
        beams;

        % frame class
        frame;

        % fem-system class
        fem;

        %%% BOUNDARY CONDITIONS FOR THE NODES %%%
        
        % list of clamped (fixed) nodes
        bc_clamped;

        % list of simply supported nodes
        bc_spprted;

        % list of infinite element nodes for harmonic analysis
        bc_infinit;

        %%% MONITORING AND PROBE NODES %%%

        % probe node indices
        probes;

        %%% CONFIGURATION %%%

        % bash (quiet) mode
        cfg_bashm = c_prj_dframe.D_CFG_BASHM;

    end

    % METHODS
    methods

        % CONSTRUCTOR
        function s = c_prj_dframe()
            
            % status message
            s.txtOut('> initializing double-frame project\n');

            s.setDefaults();
            s.initNodes();
            s.initBeams();

        end

        % SET DEFAULT VALUES
        function s = setDefaults(s)

           s.lF = s.D_LF;
           s.lInf = s.D_LINF;
           s.h = s.D_H;
           s.w = s.D_W;
           s.nSkin = s.D_NSKIN;
           s.rhoASkin = s.D_RHOASKIN;
           s.EASkin = s.D_EASKIN;
           s.EISkin = s.D_EISKIN;
           s.nFrame = s.D_NFRAME;
           s.rhoAFrame = s.D_RHOAFRAME;
           s.EAFrame = s.D_EAFRAME;
           s.EIFrame = s.D_EIFRAME;
           s.nBolt = s.D_NBOLT;
           s.rhoABolt = s.D_RHOABOLT;
           s.EABolt = s.D_EABOLT;
           s.EIBolt = s.D_EIBOLT;
           s.omega = s.D_OMEGA;
           s.f = s.D_F;
           s.fNodes = s.D_FNODES;
           s.bc_clamped = s.D_BC_CLAMPED;
           s.bc_spprted = s.D_BC_SPPRTED;
           s.bc_infinit = s.D_BC_INFINIT;
           s.probes = s.D_PROBES;

        end

        % INIALIZE FRAME CLASS
        function s = initFrame(s)

            % initialize frame-definition class
            s.frame = c_frame_def(s.nodes);

            % add the beams to frame-class
            s.frame.addBeam(s.beams);

        end

        % INITIALIZE FE-ANALYSIS
        function s = initFEA(s)

            s.initFrame();
            s.initBCs();
            s.initForces();
            s.fem = s.frame.discretize();

        end

        % INITIALIZE STRUCTURAL NODES
        function s = initNodes(s)

            % main node coordinates
            xn1 = -s.lInf;
            zn1 = 0;

            xn2 = 0;
            zn2 = zn1;

            xn3 = s.lF;
            zn3 = zn2;

            xn4 = xn3;
            zn4 = s.h;

            xn5 = xn4 + s.w;
            zn5 = zn4;

            xn6 = xn5;
            zn6 = zn3;

            xn7 = xn6 + s.lInf;
            zn7 = zn6;

            % nodes coordinate vector
            s.nodes = [ xn1 zn1 ;
                        xn2 zn2 ;
                        xn3 zn3 ;
                        xn4 zn4 ;
                        xn5 zn5 ;
                        xn6 zn6 ;
                        xn7 zn7 ];

        end

        % INITIALIZE STRUCTURAL BEAMS
        function s = initBeams(s)

            % beams specification vector
            s.beams = [ 1 2 s.rhoASkin  s.EASkin  s.EISkin  s.nSkin  ;
                        2 3 s.rhoASkin  s.EASkin  s.EISkin  s.nSkin  ;
                        3 4 s.rhoAFrame s.EAFrame s.EIFrame s.nFrame ;
                        4 5 s.rhoABolt  s.EABolt  s.EIBolt  s.nBolt  ;
                        5 6 s.rhoAFrame s.EAFrame s.EIFrame s.nFrame ;
                        6 7 s.rhoASkin  s.EASkin  s.EISkin  s.nSkin ];

        end

        % INITIALIZE BOUNDARY CONDITIONS
        function s = initBCs(s)

            % initialize clamped boundaries
            idx = s.getIndices(s.bc_clamped);
            s.frame.nodeBC_clamped(idx);

            % initialize simply supported boundaries
            idx = s.getIndices(s.bc_spprted);
            s.frame.nodeBC_jointed(idx);

            % initialize infinite element boundaries
            idx = s.getIndices(s.bc_infinit);
            s.frame.nodeBC_infinite(idx, s.omega);

        end

        % INITIALIZE NODAL FORCES
        function s = initForces(s)

            idx = s.getIndices(s.fNodes);
            s.frame.addHarmonicForce(idx, s.f);

        end

        % CALCULATE TRANSMISSION LOSS IN DOUBLE-FRAME
        function tl = calcTL(s)

            % get complex displacement amplitudes from harmonic analysis
            u = s.fem.harmonicAnalysis(s.omega);

            % get complex internal force amplitudes
            % TODO: check if this method is correct !!!
            f = s.fem.calcInnerF(u, s.omega);

            % extract the values from the probe-nodes
            u_p = s.extractNodes(u, s.probes);
            f_p = s.extractNodes(f, s.probes);

            % extract the node displacement at the force node 
            u_f = s.extractNodes(u, s.fNodes);

            % calculate powers at probe and force nodes
            P_p = sum(real(0.5*f_p.*conj(1i*s.omega*u_p)));
            P_f = sum(real(0.5*s.f.*conj(1i*s.omega*u_f)));

            % tl equals 10*log10(transmitted power through
            % double-frame/power by the excitation force)
            tl = 10*log10(abs(P_p)/abs(P_f));

        end

        % PARAMETER STUDY OF INFINITE ELEMENT BEAM LENGTH
        %
        % Inputs:
        %   p_lInf - value range for parameter lInf
        function result = study_lInf(s, p_lInf)

            % pre-allocate results vector
            result = zeros(numel(p_lInf),1);
            
            % reset project to default values
            s.setDefaults();

            % iteration counter
            c = 0;

            % begin iteration
            for lInf = p_lInf

                s.lInf = lInf;

                c = c + 1;

                % re-initialize everything
                s.initNodes();
                s.initBeams();
                s.initFEA();

                result(c) = s.calcTL();

            end % iteration

        end

        % PARAMETER STUDY OF BOLT BENDING RIGIDITY
        %
        % Inputs:
        %   p_EIBolt - value range for parameter EIBolt
        function result = study_EIBolt(s, p_EIBolt)

            % pre-allocate results vector
            result = zeros(numel(p_EIBolt),1);
            
            % reset project to default values
            s.setDefaults();

            % initialize nodes
            s.initNodes();

            % iteration counter
            c = 0;

            % begin iteration
            for EIBolt = p_EIBolt

                s.EIBolt = EIBolt;

                c = c + 1;

                % re-initialize everything
                s.initBeams();
                s.initFEA();

                result(c) = s.calcTL();

            end % iteration

        end

        % RETURN NODE INDICES OF SPECIFIC NODES (e.g. 'all')
        %
        % Inputs:
        %   p_nodeSpec - node references
        %                there are multiple options:
        %                [1 2 5]        - just returns this list of node
        %                                 indices
        %                'all'          - returns a list of ALL node
        %                                 indices
        %                'none'         - returns an empty list
        %                'firstAndLast' - returns the first and last
        %                                 node index
        %
        function idx = getIndices(s, p_nodeSpec)

            if isstr(p_nodeSpec)

                switch p_nodeSpec
                    case 'all'
                        idx = [1:size(s.nodes,1)];
                    case 'none'
                        idx = [ ];
                    case 'firstAndLast'
                        idx = [1 size(s.nodes,1)];
                    otherwise
                        idx = p_nodeSpec;
                end

            else

                idx = p_nodeSpec;

            end

        end

        % EXTRACT SPECIFIC NODE VARIABLE DATA FROM THE DISCRETIZED SYSTEM
        %
        % Inputs:
        %   p_v     - column vector with the variable data for each node
        %             of the discretized system
        %             (may be a multiple of the system node count, e.g.
        %             to extract nodal displacements or nod coordinates)
        %   p_nodes - row vector with the GLOBAL node indices which the
        %             data is extracted from
        function v = extractNodes(s, p_v, p_nodes)

            p_nodes = s.getIndices(p_nodes);

            % number of values per node
            nV = numel(p_v)/s.fem.N;

            % nodes indices of the discretized system for the main
            % system nodes (which are used in this class)
            nIdx = s.frame.frame_nodes_idx;

            % return the node values in matrix form (rows -> nodes;
            % columns -> data values)
            v = p_v(s.vectIndices(nIdx(p_nodes), nV));

        end

        % RETURN DATA VECTOR INDICES FOR SPECIFIED NODES
        %
        % Inputs:
        %   p_nIdx - node indices
        %   p_nV   - number of data values per node
        function idx = vectIndices(s, p_nIdx, p_nV)

            idx = p_nV*p_nIdx*ones(1,p_nV)-p_nV+1+ ...
                  ones(numel(p_nIdx),1)*[0:p_nV-1];

        end 

        % TEXT OUTPUT
        %
        % Inputs:
        %   p_t - text string
        function s = txtOut(s, p_t)

            if ~s.cfg_bashm
                fprintf(p_t);
            end

        end

    end

end
