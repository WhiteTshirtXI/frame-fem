function fi = beam_iForces( u, le, alpha, beam_spec )
%BEAM_IFORCES - Inner forces of a beam element @ the element nodes
% This function returns the inner forces of a beam element (normal
% force, shear force, bending moment) in global coordinates at the two
% element nodes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beam_iForces.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : fi = beam_iForces( u, le, alpha, beam_spec )
%
% Inputs :
%    u         - node displacement vector
%    le        - element length
%    alpha     - element rotation angle (rad)
%    beam_spec - beam specification containing:
%                 -> rhoA
%                 -> EA
%                 -> EI
%
% Outputs :
%    fi - inner force vector for both nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-06-09 13:25 CEST
% Last Modified : 2012-06-09 14:52 CEST
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

EIrl = EI/le;                     
EIrll = EIrl/le;

% stress-displacement matrix
%Se = @(xi) ...
%     (1/le)*[ -EA  EA          0          0           0           0 ;
%                0   0  -12*EIrll   12*EIrll      6*EIrl      6*EIrl ;
%                0   0 -6*xi*EIrl  6*xi*EIrl (3*xi-1)*EI (3*xi+1)*EI ];  
Se = (1/le) * [ -EA  EA         0         0       0       0 ;
                -EA  EA         0         0       0       0 ;
                  0   0 -12*EIrll  12*EIrll  6*EIrl  6*EIrl ;
                  0   0 -12*EIrll  12*EIrll  6*EIrl  6*EIrl ;
                  0   0    6*EIrl   -6*EIrl   -4*EI   -2*EI ;
                  0   0   -6*EIrl    6*EIrl    2*EI    4*EI ];

% calculate forces at the element nodes (xi = -1 and +1)
fi = transp(Te)*Se*Te*u;

end

