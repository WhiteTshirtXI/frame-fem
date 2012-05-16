%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frame-fem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:   vibrational analysis of a 2D-frame construction using
%                finite elements.
%
% Author:        Felix Langfeldt (felix.langfeldt@haw-hamburg.de)
%
% Creation date: 2012-05-14
% Modified:      2012-05-16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% general options %%%

% number of eigenmodes to be analysed after the elementation is done
MODAL_ANALYSIS = 4;

%%% geometric parameters %%%

% skin beam length (m)
L = 0.95;
%L = 1.00;
% double frame gap width (m)
B = 0.1;
%B = 0.0;
% frame height (m)
H = 0.1;


%%% material parameters %%%

% skin beam properties
% line mass (kg/m)
RHOA_SKIN = 27;
% axial rigidity (N)
EA_SKIN = 700e6;
% bending rigidity (N*m^2)
EI_SKIN = 583e3;

% frame properties
% line mass (kg/m)
RHOA_FRAME = 27;
% axial rigidity (N)
EA_FRAME = 700e6;
% bending rigidity (N*m^2)
EI_FRAME = 583e3;

% frame bolt properties
% line mass (kg/m)
RHOA_BOLT = 27;
% axial rigidity (N)
EA_BOLT = 700e6;
% bending rigidity (N*m^2)
EI_BOLT = 583e3;


%%% modeling parameters %%%

% number of (uniformly spaced) finite elements per structural beam element
% !! this is critical for the determination of vibrational modes etc. !!
% -> the more finite elements, the better

% skin beam
NEL_SKIN = 10;
% frame beam
NEL_FRAME = 2;
% frame bolt beam
NEL_BOLT = 2;

% minimum element length (m)
% if set to nonzero values, a check is performed, if the specified element
% numbers satisfy the constrainment of minimum element length. if not the
% element cound is altered accordingly.
LEL_MIN = 0;


% excitation frequency (rad/s)
OM = 100.0*2*pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% !!! INSERT MINIMUM LENGTH CHECK !!!


%%% MAIN BEAMS %%%

% element numbers for each beam
nel1 = NEL_SKIN;
nel2 = NEL_FRAME;
nel3 = NEL_BOLT;
nel4 = nel2;
nel5 = nel1;
%nel5 = 0;

% material properties for each beam
rhoA1 = RHOA_SKIN;
rhoA2 = RHOA_FRAME;
rhoA3 = RHOA_BOLT;
rhoA4 = rhoA2;
rhoA5 = rhoA1;

EA1 = EA_SKIN;
EA2 = EA_FRAME;
EA3 = EA_BOLT;
EA4 = EA2;
EA5 = EA1;

EI1 = EI_SKIN;
EI2 = EI_FRAME;
EI3 = EI_BOLT;
EI4 = EI2;
EI5 = EI1;

           
% beam properties vector
main_beams = [ nel1 rhoA1 EA1 EI1 ;
               nel2 rhoA2 EA2 EI2 ;
               nel3 rhoA3 EA3 EI3 ;
               nel4 rhoA4 EA4 EI4 ;
               nel5 rhoA5 EA5 EI5 ];


%%% MAIN NODES %%%

% main node coordinates
xn1 = 0;
zn1 = 0;

xn2 = L;
zn2 = zn1;

xn3 = xn2;
zn3 = H;

xn4 = xn2 + B;
zn4 = zn3;

xn5 = xn4;
zn5 = zn2;

xn6 = xn5 + L;
zn6 = zn5;

% calculations of main node indices
idx_node_1 = 1;                 % node 1
idx_node_2 = idx_node_1+nel1;   % node 2
idx_node_3 = idx_node_2+nel2;   % node 3
idx_node_4 = idx_node_3+nel3;   % node 4
idx_node_5 = idx_node_4+nel4;   % node 5
idx_node_6 = idx_node_5+nel5;   % node 6

% main node vector
main_nodes = [ xn1 zn1 idx_node_1 ;
               xn2 zn2 idx_node_2 ;
               xn3 zn3 idx_node_3 ;
               xn4 zn4 idx_node_4 ;
               xn5 zn5 idx_node_5 ;
               xn6 zn6 idx_node_6 ];

% number of nodes equals index of node 6 !
n_nodes = idx_node_6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CONSTRUCT ELEMENTS W/ ELEMENT MATRICES AND INDEX VECTORS %%%

% pre-allocate elements cell array
elements = cell(idx_node_6-1, 3);

% pre-allocate node matrix
% row structure: [n_x n_z]
nodes = zeros(n_nodes, 2);

% iteration counter
i_beam = 1;

% iterate through all main beams
for beam = main_beams'
    
    % start and end point coordinates
    xz0 = main_nodes(i_beam,1:2);
    xz1 = main_nodes(i_beam+1,1:2);
    
    % node index offset
    nOffset = main_nodes(i_beam,3);

    % min and max affected global element indices
    eMin = nOffset;
    eMax = nOffset+beam(1)-1;

    % create elements and nodes
    [elements(eMin:eMax,:), ...
     nodes(eMin:eMax+1,:)]     = create_beam(xz0, xz1, nOffset, beam);

    i_beam = i_beam + 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% ASSEMBLE SYSTEM MATRICES %%%

[MSys, KSys] = build_sys_matrix(elements);


%%% CALCULATE EIGENVALUES AND EIGENVECTORS %%%

if MODAL_ANALYSIS > 0

    % initialize figure
    fig_modes = figure;

    %number of subplots in x- and y-direction
    n_subplots_x = ceil(sqrt(MODAL_ANALYSIS));
    n_subplots_y = ceil(MODAL_ANALYSIS/n_subplots_x);

    [eigVec,eigVal] = eig(KSys, -MSys);

    % get eigenfrequencies
    eigOm = imag(sqrt(diag(eigVal)));   % rad/s
    eigF  = eigOm ./ (2*pi);            % Hz

    % mode nr.
    mode = 0;

    for i = 1:MODAL_ANALYSIS

        % if number of modes to be analyzed is greater than the actual
        % number of modes, skip for loop
        if i > size(eigF)
            break;
        end
        
        %mode = mode + 1;
        %eigF_m = eigF(end-mode-1);

        % if the eigenfrequency of the current mode is = 0, skip
        while 1
            mode = mode + 1;
            eigF_m = eigF(end-mode-1);
            if eigF_m > eps
                break
            end
        end

        % get mode shape
        nodes_m = get_mode_shape(nodes, eigVec(:,end-mode-1));

        % plot mode
        fig_modes_sub = subplot(n_subplots_y, n_subplots_x, i, ...
                                'parent', fig_modes);
        plot_mode_shape(nodes, nodes_m, i, eigF_m);

        % info output
        fprintf('Mode %i : f = %8.2f Hz\n', [i, eigF_m]);


    end

end


%%% TRANSIENT ANALYSIS %%%

% assemble matrices for transient ode-system
MSys_y = [             MSys zeros(3*n_nodes) ;
           zeros(3*n_nodes)   eye(3*n_nodes) ];
MSys_y = sparse(MSys_y);

KSys_y = [ zeros(3*n_nodes)            KSys ;
            -eye(3*n_nodes) zeros(3*n_nodes) ];
KSys_y = sparse(KSys_y);


% create force amplitude vector
f_a = zeros(3*n_nodes,1);
% the first node gets a z-force amplitude
f_a(2) = 100.0;
% expand force amplitude vector
f_a = [f_a ; zeros(3*n_nodes,1)];

% force function vector
f_y = @(t) f_a*real(exp(j*OM*t));

% right hand side of system ODE
f_rhs = @(t,y) f_y(t) - KSys_y*y;

% initial values (all zero!)
y_0 = zeros(2*3*n_nodes,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot_nodes(nodes);
