% [u0, x, y, N] = generate2DHexLattice(S,alpha)
%
% Generates an initial guess (u0) for equilibrium positions for the number
% of shells given by shells (based on triangular/hexagonal lattice structure)
%
% INPUT
% S:      Number of complete hex shells to add
% alpha:  Dimensionless lattice spacing (lattice constant)
% 
% OUTPUT
% x, y:   Lattice positions in each direction
% u0:     Single array of these positions with format
%         [x1 x2 ... xN ... y1 y2 ... yN]
% N:      total number of ions for S shells

function [u0, x, y, N] = generate2DHexLattice(S,alpha)

u0 = [0 0];    		% Center point always 0,0

for s=1:S     		% Add next shell
    u0 = [u0 addHexShell(s)];  
end
u0 = u0*alpha; 		% Set lattice constant (addHexShell has lattice constant = 1)
N = 1 + 6*sum(1:S); % Number of ions in configuration
x = u0(1:2:end);    % Format: [x1 y1 x2 y2 ... xN yN]
y = u0(2:2:end);    % It was only in this format to add multiple shells

%REFORMAT to [x1 x2 ... xN ... y1 y2 ... yN]
u0(1:N) = x;   
u0(N+1:end) = y;

end

        
   