function [ MSys, KSys ] = build_sys_matrix( elData )
%BUILD_SYS_MATRIX - Assemble system matrices using element data
% This function returns the system mass and stiffness matrices using the
% element data provided as a parameter. The element data vector contains
% the element matrices and the index vectors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build_sys_matrix.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Syntax : [ MSys, KSys ] = build_sys_matrix( elData )
%
% Inputs :
%    elData - A cell matrix containing in each row:
%              -> element mass matrix
%              -> element stiffness matrix
%              -> element index vector
%
% Outputs :
%    MSys - System mass matrix (sparse)
%    KSys - System stiffness matrix (sparse)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author        : Felix Langfeldt
%                 felix.langfeldt@haw-hamburg.de
%
% Creation Date : 2012-05-17 12:28 CEST
% Last Modified : 2012-05-18 11:20 CEST
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if the cell arrays of element matrices and element index vectors
% are equally sized

Me = elData(:,1);
Ke = elData(:,2);
vecIe = elData(:,3);

if(numel(Me) == numel(vecIe) && numel(Ke) == numel(vecIe))

    % calculate system matrix bandwidth
    % to do that, we have to find out the maximum index difference of
    % the index vectors
    maxDiff = 0;
    for indexVector = vecIe'
        delta =	max(indexVector{:}) - min(indexVector{:});
        maxDiff = max(maxDiff, delta);
    end

    % max. bandwith equals maxDiff + 1
    bandWith = maxDiff + 1;  
    
    % find out maximum matrix index -> number of nodes
    N = max(cat(2,vecIe{:}));
    
    MSys = spalloc(N,N,(2*bandWith+1)*N);
    KSys = spalloc(N,N,(2*bandWith+1)*N);
    
    % run through all elements
    for e = 1:numel(Me)
        
        % element matrices
        Me_e = Me{e};
        Ke_e = Ke{e};
        % index vector
        vecIe_e = vecIe{e};
        
        % run through element matrix elements
        for i = 1:size(Me_e,1)
            for j = 1:size(Me_e,2)
                MSys(vecIe_e(i),vecIe_e(j)) = MSys(vecIe_e(i),vecIe_e(j)) + Me_e(i,j);
                KSys(vecIe_e(i),vecIe_e(j)) = KSys(vecIe_e(i),vecIe_e(j)) + Ke_e(i,j);
            end
        end
        
    end
        
end

end

