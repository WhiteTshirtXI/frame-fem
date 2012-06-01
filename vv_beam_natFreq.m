function [fn] = vv_beam_natFreq( p_n, p_nSimply, p_nClamped )
%VV_BEAM_NATFREQ - Natural frequencies of a simple beam
% This functions returns the natural frequencies of a simple, homogenous
% euler-beam of unit length, unit bending rigidity EI and unit mass ber
% length rhoA for a specified number of modes and specified boundary
% conditions. The reference used is:
%
% YOUNG, D. and R. FELGAR: Tables of characteristic functions
%                          representing normal modes of vibration of a
%                          beam. Engineering Research Series. University
%                          of Texas, 1949.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vv_beam_natFreq.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : [fn] = vv_beam_natFreq( p_n, p_nSimply, p_nClamped )
%
% Inputs :
%    p_n        - number of modes
%    p_nSimply  - number of simply supported nodes (0, 1 or 2)
%    p_nClamped - number of clamped nodes (0, 1 or 2)
%
% Outputs :
%    fn - vector of p_n eigenfrequencies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-31 10:37 CEST
% Last Modified : 2012-05-31 15:19 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pre-allocate frequency vector
    fn = zeros(p_n,1);

    % in the reference there are numerical values for the first five
    % natural frequencies tabulated. for higher frequencies, analytical
    % solutions are given.

    % first check if the overall number of bc-specifications doesn't
    % exceed 2, because the beam only has two nodes
    if sum([p_nSimply p_nClamped]) <= 2

        switch p_nSimply

            % CASE 1 : supported-supported beam
            case 2

                fn_5 = [  1*pi ;
                          2*pi ;
                          3*pi ;
                          4*pi ;
                          5*pi ].^2;

                % inline function for frequencies > mode 5
                fn_f = inline('(n*pi).^2');

            % CASE 2 : clamped-supported or free-supported beam
            case 1

                fn_5 = [  3.92660230 ;
                          7.06858275 ;
                         10.21017613 ;
                         13.35176878 ;
                         16.49336143 ].^2;

                % inline function for frequencies > mode 5
                fn_f = inline('((4*n+1)*pi/4).^2');

            % no simply supported bc
            case 0

                switch p_nClamped

                    % CASE 3 : clamped-free beam
                    case 1

                        fn_5 = [  1.87510410 ;
                                  4.69409113 ;
                                  7.85475743 ;
                                 10.99554074 ;
                                 14.13716839 ].^2;

                        % inline function for frequencies > mode 5
                        fn_f = inline('((2*n-1)*pi/2).^2');

                    % CASE 4 : clamped-clamped or free-free beam
                    otherwise

                        fn_5 = [  4.73004080 ;
                                  7.85320460 ;
                                 10.99560730 ;
                                 14.13716550 ;
                                 17.27875960 ].^2;

                        % inline function for frequencies > mode 5
                        fn_f = inline('((2*n+1)*pi/2).^2');

                end

        end

        % get first max. 5 values from tabulated values
        nFreq_5 = min(5,p_n);

        fn(1:nFreq_5) = fn_5(1:nFreq_5);

        % get the rest (if neccessary) from the inline functions
        nFreq_rest = [1:p_n-nFreq_5]'+5;

        fn(nFreq_rest:end) = fn_f(nFreq_rest);

    end

end
