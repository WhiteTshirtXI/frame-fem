%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vv_beam_inf.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% =================================================
% | V&V - BEAM (FORCED VIBRATIONS, SEMI-INFINITE) |
% =================================================
%
% Description   : verification of the fem-system class using the forced
%                 harmonic vibrations of a semi-infinite beam which has
%                 exact analytical solutions as a reference. Beam
%                 properties are taken from:
%
%                 WANG, C. and J.C.S. LAI (2000): Modelling the
%                 Vibration Behaviour of Infinite Structures By FEM.
%                 Journal of Sound and Vibration 229(3), pp. 453-466.
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-06-19 14:40 CEST
% Last Modified : 2012-06-19 14:48 CEST
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

% angular frequencies for harmonic analysis
OM = [100]*2*pi;

% excitation force amplitude (N)
F = 100.0;

%%% beam properties (materials, geometry, etc.) %%%

% square crossection - height and width (m)
H = 0.001;
W = 0.01;
% density (kg/m^3)
RHO = 7860;
% young's modulus (N/m^2)
E = 210000e6;

% length of beam (m)
L = 0.2;
% beam rotation angle (rad)
ALPHA = 0*pi/180;

% number of finite elements
% (use multiple values for parameter study)
NEL = [20];


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


% bending wave velocities for all frequencies
cB = (EI/RHOA)^.25*sqrt(OM);


%%% FRAME STRUCTURE %%%

%%% beam nodes %%%
beam_nodes = L*[ 0          0          ;
                 cos(ALPHA) sin(ALPHA) ];


                 
%%% HARMONIC ANALYSIS OF SEMI-INFINITE BEAM %%%
fprintf('HARMONIC ANALYSIS OF SEMI-INFINITE BEAM\n');

%%% do the following for every entry in the element number vector %%%

% run number
i_run = 0;

for i_nel = NEL

    i_run = i_run+1;

    fprintf('RUN #%02i : %2i element(s)\n', [i_run,i_nel]);

    %%% run through all excitation frequencies %%%
    for i_om = [1:numel(OM)]

        o = OM(i_om);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% create frame class %%%
        frame = c_frame_def(beam_nodes);

        % add beam to frame
        frame.addBeam( [ 1 2 RHOA EA EI i_nel/L]);

        % add harmonic force to the first node
        frame.addHarmonicForceZ(1, F);

        % apply infinite boundary condition at the end node
        frame.nodeBC_infinite(2, o);

        % discretize the system
        sys_fem = frame.discretize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% HARMONIC ANALYSIS %%%

        fprintf(['  f = %8.2f Hz - '], o/(2*pi));
       
        % calculate vector of complex displacement amplitudes for the
        % given frequency
        u = sys_fem.harmonicAnalysis(o);

        % get displacement amplitude at the first node in z-direction
        uz1 = u(3*(1-1)+2);

        % calculate velocity amplitude
        vz1 = 1i*o*uz1;

        % real part of the driving point mobility at the excitation
        % point (node 1)
        dpm1_calc = full(real(vz1/F));

        % reference driving point mobility according to CREMER and HECKL
        % (1967)
        dpm1_ref = real(1/(0.5*RHOA*cB(i_om)*(1+1i)));

        % relative difference
        dpm1_diff = (dpm1_calc-dpm1_ref)/dpm1_ref;

        fprintf(['driving point mobilities :\n' ...
                 '                    ... calculated : %6.2e m/(Ns)\n'...
                 '                    ... reference  : %6.2e m/(Ns)\n'...
                 '                    ... difference : %8.3f %%\n' ...
            ], [dpm1_calc, dpm1_ref, dpm1_diff*100]);


        fprintf('\n');

    end

end

