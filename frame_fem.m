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
% Last Modified : 2012-06-06 16:06 CEST
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
MODAL_ANALYSIS = 4;

% activate plot (show = 1, hide = 0)
show_plots = 1;


% skin beam length (m)
L = 3;
%L = 1.00;
% double frame gap width (m)
B = 0.1;
%B = 0.0;
% frame height (m)
H = 0.12;


%%% material parameters %%%

% skin beam properties
% line mass (kg/m)
RHOA_SKIN = 5.54;
% axial rigidity (N)
EA_SKIN = 143.7e6;
% bending rigidity (N*m^2)
EI_SKIN = 12266.5;

% frame properties
% line mass (kg/m)
RHOA_FRAME = 1.17;
% axial rigidity (N)
EA_FRAME = 30.24e6;
% bending rigidity (N*m^2)
EI_FRAME = 66949.1 ;

% frame bolt properties
% line mass (kg/m)
RHOA_BOLT = 0.788;
% axial rigidity (N)
EA_BOLT = 20.5e6;
% bending rigidity (N*m^2)
EI_BOLT = 247.397;


%%% modeling parameters %%%

% number of (uniformly spaced) finite elements per structural beam element
% !! this is critical for the determination of vibrational modes etc. !!
% -> the more finite elements, the better

% skin beam
NEL_SKIN = 30;
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

% status message
fprintf('> pre-processing ...\n');

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

% frame nodes coordinate vector
frame_nodes = [ -L   0  ;
                xn1 zn1 ;
                xn2 zn2 ;
                xn3 zn3 ;
                xn4 zn4 ;
                xn5 zn5 ;
                xn6 zn6 ];

% frame beams specification vector
frame_beams = [ 1 2 rhoA1 EA1 EI1 nel1 ;
                2 3 rhoA1 EA1 EI1 nel1 ;
                3 4 rhoA2 EA2 EI2 nel2 ;
                4 5 rhoA3 EA3 EI3 nel3 ;
                5 6 rhoA4 EA4 EI4 nel4 ;
                6 7 rhoA5 EA5 EI5 nel5 ];


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


% status message
fprintf('  ... discretized the system with %i nodes\n', sys_fem.N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% status message
fprintf('\n> processing ...\n');


%%% ASSEMBLE SYSTEM MATRICES %%%

MSys = sys_fem.MSys;
KSys = sys_fem.KSys;

%%% CALCULATE EIGENVALUES AND EIGENVECTORS %%%

if MODAL_ANALYSIS > 0

    % status message
    fprintf(['  ... calculating the %i lowest normal modes of the'   ...
             ' system\n'], MODAL_ANALYSIS);

    if show_plots == 1
        sys_fem.plotModes(MODAL_ANALYSIS);

        try
            waitforbuttonpress();
            close;
        end
        
    end
end



% status message
fprintf(['  ... performing a harmonic analysis for the system'   ...
         ' @ f = %6.1f Hz\n'], OM/(2*pi));

%%% HARMONIC ANALYSIS %%%

% add harmonic force of 100 N in z-direction to the first node
frame.addHarmonicForceZ(2, 100.0);

% add infinite boundary element to the end nodes
frame.nodeBC_infinite([1 7], OM);

% discretize the system
sys_fem = frame.discretize();

% perform harmonic analysis to calculate the vector of complex
% displacement amplitudes uH_c
uH_c = sys_fem.harmonicAnalysis(OM);


% POWER FLOW ANALYSIS
% calculate complex displacement velocity amplitudes
vH_c = 1i*OM*uH_c;

% calculate complex internal forces at each node
fi_c = sys_fem.calcInnerF(uH_c, OM);

% calculate power at each node
Pn = zeros(sys_fem.N, 1);
for n = [1:sys_fem.N]

    % index of first dof at current node
    n0 = sys_fem.NDOF*(n-1)+1;
    % index of last dof at current node
    n1 = n0+sys_fem.NDOF-1;

    % sum over all dofs
    Pn(n) = 0.5*real(vH_c(n0:n1)'*fi_c(n0:n1));

end

% rearrange the complex displacement amplitude vector to yield the x-
% and z-components
uH_c_xz = [uH_c(1:3:end) uH_c(2:3:end)];

% get real part of the complex displacement amplitudes
uH_xz = real(uH_c_xz);

% get magnitude of the complex displacement amplitudes
%uH_xz = abs(uH_c_xz);

% normalize for nicer plotting
uH_xz = 0.1*uH_xz/max(abs(uH_xz(:)));

% visualize displacement amplitudes
if show_plots == 1;
    sys_fem.plot_nodes.plotDisplaced(sys_fem.nodes, uH_xz, sys_fem.nAdj, ...
            {sprintf('Harmonic Analysis: f = %8.2f Hz', OM/(2*pi))});
    try
        waitforbuttonpress();
        close;
    end
end

% calculate Transmission loss (dB) by the ratio of power input and power
% output
TL=10*log10((abs(Pn(end)))/(abs(Pn(1))+abs(Pn(end))));
% status message of Transmission loss 
fprintf('\n%4.1f dB\n',full(TL));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

