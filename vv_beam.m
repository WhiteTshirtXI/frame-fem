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
% Last Modified : 2012-05-31 16:13 CEST
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

% BENDING MODES ONLY -> axial rigidity to zero
EA = 0;


%%% beam nodes %%%
beam_nodes = L*[ 0          0          ;
                 cos(ALPHA) sin(ALPHA) ];

%%% BOUNDARY CONDITIONS %%%
% clamped nodes 
nodes_clamped = [ 1 ];
% jointed nodes
nodes_jointed = [ 2 ];


%%% according to the lowest element number and the selected boundary %%%
%%% conditions, the maximum number of modes might need a correction  %%%

% corrected maximum modes count = (minimum number of elements + 1) *
%                                 2 beam node DOFs - constrained dofs
%                                 because of BCs
maxModes_corrected = (min(NEL)+1)*2-numel(nodes_clamped)*2           ...
                                   -numel(nodes_jointed);

if maxModes_corrected < MAXMODES
    fprintf(['CORRECTING maximum modes number to %i ...\n\n'],       ...
             maxModes_corrected);
    MAXMODES = maxModes_corrected;
end


% get reference frequencies
f_n_ref = vv_beam_natFreq(MAXMODES, numel(nodes_jointed),            ...
                          numel(nodes_clamped))                      ...
          *sqrt(EI/RHOA)/(L^2*2*pi);


%%% create eigenfrequency matrix where the rows correspond to each   %%%
%%% mode and the columns correspond to each refinement stage         %%%
f_results = zeros(MAXMODES, numel(NEL));


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

        % plot eigenmodes
        sys_fem.plotModes(MAXMODES);

        % calculate eigenfrequencies
        eigF = sort(sys_fem.eigF(MAXMODES));

        % store eigenfrequencies in results matrix
        f_results(:,i_run) = eigF;

        % calculate difference from reference frequencies
        f_diff = (eigF-f_n_ref)./f_n_ref;

        % mode number
        m = 1;

        % write out all eigenfrequencies
        for f = eigF'

            fprintf(['  ... mode #%02i :      f = %8.1f Hz\n'   ...
                     '                    ref = %8.1f Hz\n'   ...
                     '                   diff = %8.3f %%\n'], ...
                    [m,f,f_n_ref(m),f_diff(m)]);

            m = m + 1;

        end

    end

    fprintf('\n');

end
