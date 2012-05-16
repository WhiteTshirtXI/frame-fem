function [ el, nodes ] = create_beam( xz0, xz1, offset, beam_specs )
%CREATE_BEAM Creates a beam with specified number of elements (uniformly
%spaced)
%   Input arguments:  xz0        - start point coordinate
%                     xz1        - end point coordinate
%                     offset     - node index offset
%                     beam_specs - beam properties vector containing:
%                                  [Ne rhoA EA EI]
%   Output arguments: el         - cell array containing Ne rows with:
%                                  [Me Ke idxVec]
%                     nodes      - node data matrx containing Ne+1 rows:
%                                  [n_x n_y]

% calculate beam length and angle
dxz = xz1 - xz0;
L = norm(dxz);
alpha = -atan2(dxz(2),dxz(1));  % negative because of coordinate system!

% number of elements
Ne = beam_specs(1);

% calculate element length
le = L/Ne;

% calculate element matrices
[Me, Ke] = matrices_beam(le, alpha, beam_specs(2:4));

% pre-allocate results cell array
el = cell(Ne,3);

% pre-allocate nodes matrix
nodes = zeros(Ne+1,2);
nodes(1,:) = xz0;

% create index vectors for each element
for e = 1:Ne
    % first node index
    k1 = (offset-1)+e;
    % second node index
    k2 = k1+1;
    
    % compose index vector
    k1i = 1+3*(k1-1);
    k2i = 1+3*(k2-1);
    idxVec = [k1i k2i k1i+1 k2i+1 k1i+2 k2i+2];
    
    % add to results cell array
    el(e,:) = {Me  Ke  idxVec};

    % add coordinates to nodes matrix
    nodes(e+1,:) = xz0 + e*le*[cos(-alpha) sin(-alpha)];
    
end


end

