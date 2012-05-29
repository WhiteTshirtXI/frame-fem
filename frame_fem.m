%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frame_fem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% =============
% | FRAME-FEM |
% =============
%
% Description   : vibrational analysis of a 2D-frame construction using
%                 finite elements.
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-14 14:00 CEST
% Last Modified : 2012-05-29 17:26 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% general options %%%

% number of eigenmodes to be analysed after the elementation is done
MODAL_ANALYSIS = 6;

% Materials

% Alu
% Density (kg/m3)
RHO_ALU=2700; 
% Youngs Modulus (MPa)
E_ALU=70e3; 
% poisson constant (-)
MUE_ALU=0.33;
% Steel
% Density (kg/m3)
RHO_STEEL=8100;
% Youngs Modulus (MPa)
E_STEEL=210e3; 
% poisson constant (-)
MUE_STEEL=0.3;


%%% geometric parameters %%%

% skin beam properties
% skin height (m)
SKIN_HEIGHT=0.013;
% skin width (m)
SKIN_WIDTH=0.68;     
% skin length (m)
SKIN_LENGTH=1.2;   

% frame properties
% frame height (m)
FRAME_HEIGHT=0.12;  
% frame width (m)
FRAME_WIDTH=SKIN_WIDTH;  
% frame length (m)
FRAME_LENGTH=0.004;  

% frame bolt properties
% bolt diameter (m)
BOLT_DIAMETER=0.01;
% bolt length (m) = double frame gap width (m)
BOLT_LENGTH=0.12;


%%% material parameters %%%

% skin beam properties
% line mass (kg/m), % axial rigidity (N), % bending rigidity (N*m^2)
[ rhoA_skin,EA_skin,EI_skin ] = prop_sec_rect(RHO_ALU,E_ALU,MUE_ALU,SKIN_LENGTH,SKIN_WIDTH,SKIN_HEIGHT);

% frame properties
% line mass (kg/m), % axial rigidity (N), % bending rigidity (N*m^2)
[ rhoA_frame,EA_frame,EI_frame ] = prop_sec_rect(RHO_ALU,E_ALU,MUE_ALU,FRAME_HEIGHT,FRAME_WIDTH,FRAME_LENGTH);

% frame bolt properties
% line mass (kg/m), % axial rigidity (N), % bending rigidity (N*m^2)
[ rhoA_bolt,EA_bolt,EI_bolt ] = prop_sec_circle(RHO_ALU,E_ALU,MUE_ALU,BOLT_LENGTH,BOLT_DIAMETER);

% wieder Löschen!
rhoA_bolt=rhoA_frame;
EA_bolt=EA_frame;
EI_bolt=EI_frame;


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
rhoA1 = rhoA_skin;
rhoA2 = rhoA_frame;
rhoA3 = rhoA_bolt;
rhoA4 = rhoA2;
rhoA5 = rhoA1;

EA1 = EA_skin;
EA2 = EA_frame;
EA3 = EA_bolt;
EA4 = EA2;
EA5 = EA1;

EI1 = EI_skin;
EI2 = EI_frame;
EI3 = EI_bolt;
EI4 = EI2;
EI5 = EI1;


%%% MAIN NODES %%%

% main node coordinates
xn1 = 0;
zn1 = 0;

xn2 = SKIN_LENGTH;
zn2 = zn1;

xn3 = xn2;
zn3 = FRAME_HEIGHT;

xn4 = xn2 + BOLT_LENGTH;
zn4 = zn3;

xn5 = xn4;
zn5 = zn2;

xn6 = xn5 + SKIN_LENGTH;
zn6 = zn5;

% frame nodes coordinate vector
frame_nodes = [ xn1 zn1 ;
                xn2 zn2 ;
                xn3 zn3 ;
                xn4 zn4 ;
                xn5 zn5 ;
                xn6 zn6 ];

% frame beams specification vector
frame_beams = [ 1 2 rhoA1 EA1 EI1 nel1 ;
                2 3 rhoA2 EA2 EI2 nel2 ;
                3 4 rhoA3 EA3 EI3 nel3 ;
                4 5 rhoA4 EA4 EI4 nel4 ;
                5 6 rhoA5 EA5 EI5 nel5 ];


%%% BOUNDARY CONDITIONS %%%
% clamped nodes 
nodes_clamped = [ ];
% jointed nodes
nodes_jointed = [ ];


% create frame-class
frame = c_frame_def(frame_nodes);

% add beams
frame.addBeam(frame_beams);

% apply boundary conditions
frame.nodeBC_clamped(nodes_clamped);
frame.nodeBC_jointed(nodes_jointed);

% discretize the system
sys_fem = frame.discretize();

% number of nodes 
n_nodes = sys_fem.N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% ASSEMBLE SYSTEM MATRICES %%%

MSys = sys_fem.MSys;
KSys = sys_fem.KSys;

%%% CALCULATE EIGENVALUES AND EIGENVECTORS %%%

if MODAL_ANALYSIS > 0

    sys_fem.plotModes(MODAL_ANALYSIS);

end


%%% HARMONIC ANALYSIS %%%

% add infinite boundary element to the end nodes
frame.nodeBC_infinite(6, OM);

% discretize the system
sys_fem = frame.discretize();

%%% TRANSIENT ANALYSIS %%%

% initialize sparse matrices for transient ode-system
% -> (2*number of dof)x(2*number of dof)-matrices
MSys_y = sparse(2*3*n_nodes,2*3*n_nodes);
KSys_y = sparse(2*3*n_nodes,2*3*n_nodes);


% the following steps are important for keeping the extended system
% matrices for the transient ode-system as sparse as possible !!!

% assemble extended mass matrix in two steps:
% STEP 1: for every odd row and column of MSys_y add the
%         corresponding matrix element of the system mass matrix
%         (this leads to a matrix of the following structure
%                    | M11  0  M12  0  M13  0 |
%                    |  0   0   0   0   0   0 |
%                    | M21  0  M22  0  M23  0 |
%                    |  0   0   0   0   0   0 |
%                    | M31  0  M32  0  M33  0 | 
%                    |  0   0   0   0   0   0 | )
MSys_y(1:2:end,1:2:end) = MSys; 

% STEP 2: fill the even diagonal elements with ones 
%         (this finally leads to a matrix of the following structure
%                    | M11  0  M12  0  M13  0 |
%                    |  0   1   0   0   0   0 |
%                    | M21  0  M22  0  M23  0 |
%                    |  0   0   0   1   0   0 |
%                    | M31  0  M32  0  M33  0 |
%                    |  0   0   0   0   0   1 | )
MSys_y(2:2:end,2:2:end) = speye(3*n_nodes);

% assemble extended stiffness matrix in two steps:
% STEP 1: for every odd row and even column of KSys_y add the
%         corresponding matrix element of the system stiffness matrix
%         (this leads to a matrix of the following structure
%                    | 0  K11  0  K12  0  K13 |
%                    | 0   0   0   0   0   0  |
%                    | 0  K21  0  K22  0  K23 |
%                    | 0   0   0   0   0   0  |
%                    | 0  K31  0  K32  0  K33 |
%                    | 0   0   0   0   0   0  | )
KSys_y(1:2:end,2:2:end) = KSys; 

% STEP 2: fill the even elements of the first lower diagonal with
%         negative ones
%         (this finally leads to a matrix of the following structure
%                    | 0  K11  0  K12  0  K13 |
%                    | -1  0   0   0   0   0  |
%                    | 0  K21  0  K22  0  K23 |
%                    | 0   0   -1  0   0   0  |
%                    | 0  K31  0  K32  0  K33 |
%                    | 0   0   0   0   -1  0  | )
KSys_y(2:2:end,1:2:end) = -speye(3*n_nodes);


% initialize sparse force amplitude vector
% -> (2*number of dof)-vector
f_a = sparse(2*3*n_nodes,1);
% the first node gets a z-force amplitude
% position in force vector: 1+6*(node_nr-1)+2*(dof_nr-1)
f_a(3) = 100.0;

% force function vector
f_y = @(t) f_a*imag(exp(j*OM*t));

% right hand side of system ODE
f_rhs = @(t,y) f_y(t) - KSys_y*y;

% initial values (all zero!)
y_0 = sparse(2*3*n_nodes,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

