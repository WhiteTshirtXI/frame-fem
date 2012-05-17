function [ Me, Ke ] = matrices_beam( le, alpha, beam_spec )
%MATRICES_BEAM - Return element mass and stiffnes matrices of a beam
% This function returns the element mass and stiffness matrices of a
% truss-beam-element which the frame-structures are composed of.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrices_beam.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : [ Me, Ke ] = matrices_beam( le, alpha, beam_spec )
%
% Inputs :
%    le        - element length
%    alpha     - element rotation angle (rad)
%    beam_spec - beam specification containing:
%                 -> rhoA
%                 -> EA
%                 -> EI
%
% Outputs :
%    Me - coordinate transformed element mass matrix
%    Ke - coordinate transformed element stiffness matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-17 12:39 CEST
% Last Modified : 2012-05-17 12:43 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% standard shape functions of trusses/beams according to LIU/QUEK 2003
              
ca = cos(alpha);
sa = sin(alpha);

% transformation matrix
              
Te = [ ca  0 -sa   0 0 0 ;
        0 ca   0 -sa 0 0 ;
       sa  0  ca   0 0 0 ;
        0 sa   0  ca 0 0 ;
        0  0   0   0 1 0 ;
        0  0   0   0 0 1 ];
    
% get beam specification
rhoA = beam_spec(1);
EA = beam_spec(2);
EI = beam_spec(3);


Me = (rhoA*le/420) * [ 140  70      0      0       0       0 ;
                        70 140      0      0       0       0 ;
                         0   0    156     54  -22*le   13*le ;
                         0   0     54    156  -13*le   22*le ;
                         0   0 -22*le -13*le  4*le^2 -3*le^2 ;
                         0   0  13*le  22*le -3*le^2  4*le^2 ];

EIrl = EI/le;                     
EIrll = EIrl/le;
                     
Ke = (1/le) * [  EA -EA         0         0       0       0 ;
                -EA  EA         0         0       0       0 ;
                  0   0  12*EIrll -12*EIrll -6*EIrl -6*EIrl ;
                  0   0 -12*EIrll  12*EIrll  6*EIrl  6*EIrl ;
                  0   0   -6*EIrl    6*EIrl    4*EI    2*EI ;
                  0   0   -6*EIrl    6*EIrl    2*EI    4*EI ];
              
% transform matrices
Me = transp(Te)*Me*Te;
Ke = transp(Te)*Ke*Te;


end

