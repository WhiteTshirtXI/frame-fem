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
% Last Modified : 2012-05-21 13:45 CEST
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

% square crossection - height and width (m)
H = 0.1;
W = 0.1;
% density (kg/m^3)
RHO = 2700;
% young's modulus (N/m^2)
E = 72000e6;

% length of beam (m)
L = 1;
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

    % calculate element length
    l_el = L/i_nel;

    % delta x and delta z for each node
    deltaxz_node = l_el*[cos(ALPHA) -sin(ALPHA)];

    % create nodes
    nodes = [0:i_nel]'*deltaxz_node;

    % create element node indices vector
    % (this is rather simple: just put an element between each
    % neighbored node)
    nodes_el = [[1:i_nel]' [2:i_nel+1]'];

    % initialize sys_fem-class using the node vector
    sys_fem = c_sys_fem(nodes);

    % add elements and assemble system matrices
    sys_fem.add_element(nodes_el, RHOA, EA, EI);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% MODAL ANALYSIS %%%
    
    if MAXMODES > 0

        fprintf(['  Performing a modal analysis for the max. %i ' ...
                    'lowest eigenmodes ...\n'], MAXMODES);

        % get eigenvalues and eigenvectors
        [eig_val,eig_vec] = sys_fem.getModes(MAXMODES);

        % calculate eigenfrequencies
        eig_frq = imag(sqrt(eig_val))/(2*pi);

        % write out all eigenfrequencies
        for m = 1:numel(eig_frq)

            fprintf('  ... mode #%02i : f = %8.2f Hz\n', ...
                    [m,eig_frq(m)]);

        end

    end

    fprintf('\n');

end
