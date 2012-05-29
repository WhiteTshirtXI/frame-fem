%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vv_beam.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ==============
% | V&V - BEAM |
% ==============
%
% Description   : verification of the fem-system class using a simple
%                 beam (cantilever) structure.
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-21 09:51 CEST
% Last Modified : 2012-05-29 09:55 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% CLEAN UP WORKSPACE %%%
clear all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% analysis settings %%%

% max. number of eigenmodes
MAXMODES = 4;

%%% beam properties (materials, geometry, etc.) %%%
%%% according to FAHY 2007                      %%%

% square crossection - height and width (m)
H = 0.005;
W = 0.02;
% density (kg/m^3)
RHO = 2700;
% young's modulus (N/m^2)
E = 71000e6;

% length of beam (m)
L = 0.5;
% beam rotation angle (rad)
ALPHA = 0*pi/180;

% number of finite elements
% (use multiple values for parameter study)
NEL = [1 2];


%%% derived properties %%%

% crossectional area (m^2)
A = H*W;
% area moment of inertia (m^4)
I = W*H^3/12;
% line mass (kg/m)
RHOA = RHO*A;
% axial rigidity (N)
EA = E*A;
% bending rigidity (N*m^2)
EI = E*I;


%%% beam nodes %%%
beam_nodes = L*[ 0          0          ;
                 cos(ALPHA) sin(ALPHA) ];

%%% BOUNDARY CONDITIONS %%%
% clamped nodes 
nodes_clamped = [ 1 ];
% jointed nodes
nodes_jointed = [ 2 ];


%%% do the following for every entry in the element number vector %%%

% run number
i_run = 0;

for i_nel = NEL

    i_run = i_run+1;

    fprintf('RUN #%02i : %2i element(s)\n', [i_run,i_nel]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% create frame class %%%
    frame = c_frame_def(beam_nodes);

    % add beam to frame
    frame.addBeam( [ 1 2 RHOA EA EI i_nel]);

    % apply boundary conditions
    frame.nodeBC_clamped(nodes_clamped);
    frame.nodeBC_jointed(nodes_jointed);

    % discretize the system
    sys_fem = frame.discretize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% MODAL ANALYSIS %%%
    
    if MAXMODES > 0

        fprintf(['  Performing a modal analysis for the max. %i ' ...
                    'lowest eigenmodes ...\n'], MAXMODES);

        sys_fem.plotModes(MAXMODES);


%        % write out all eigenfrequencies
%        for m = 1:numel(eig_frq)
%
%            fprintf('  ... mode #%02i : f = %8.2f Hz\n', ...
%                    [m,eig_frq(m)]);
%
%        end

    end

    fprintf('\n');

end
