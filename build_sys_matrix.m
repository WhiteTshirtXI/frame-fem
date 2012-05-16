function [ Msys Ksys ] = build_sys_matrix( elData )
%BUILD_SYS_MATRIX Build system matrix from element matrices and index
%vectors
%   Input parameters:  elData - a cell matrix containing in each row:
%                               -> element mass matrix
%                               -> element stiffness matrix
%                               -> element index vector
%   Output parameters: Msys   - system matrix

% check if the cell arrays of element matrices and element index vectors
% are equally sized

Me = elData(:,1);
Ke = elData(:,2);
vecIe = elData(:,3);

if(numel(Me) == numel(vecIe) && numel(Ke) == numel(vecIe))
    
    % find out maximum matrix index -> number of nodes
    N = max(cat(2,vecIe{:}));
    
    Msys = zeros(N);
    Ksys = Msys;
    
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
                Msys(vecIe_e(i),vecIe_e(j)) = Msys(vecIe_e(i),vecIe_e(j)) + Me_e(i,j);
                Ksys(vecIe_e(i),vecIe_e(j)) = Ksys(vecIe_e(i),vecIe_e(j)) + Ke_e(i,j);
            end
        end
        
    end
        
end

end

