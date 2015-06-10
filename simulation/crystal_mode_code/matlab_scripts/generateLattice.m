% [u0, x, y, N] = generateLattice(Nf,alpha)
%
% Generates an initial guess lattice (u0) to solve for equilibrium positions for Nf ions
%
% This function finds the maximum integer number of complete subshells that can be generated
% with Nf ions and generates that lattice as the base lattice (generate2DHexLattice.m).
% It then generates the next possible complete subshell.
% For each of the remaining ions to add, the potential energy of adding each ion from the
% next possible complete subshell is calculated, and the ion with the smallest increase in 
% potential energy is added to the crystal permanently
%
% INPUT
% Nf:     Number of ions in final configuration
%		  NOTE: Unless you would like to study a particular value of N
%             I recommend choosing N = 1 + 6*sum(1:S); for S an integer
%			  as this code is not very efficient and slow for large N not of the form above
%
% alpha:  Dimensionless lattice spacing (lattice constant)
% 		  (If you don't know the approximate scaling, set alpha = 1.
%          After this funciton generates a lattice, call idealLatticeScale(...) )
% 
% OUTPUT
% x, y:   Lattice positions in each direction
% u0:     Single array of these positions with format
%         [x1 x2 ... xN ... y1 y2 ... yN]

function [u0, x, y] = generateLattice(Nf,alpha)

global N M 
% Nf is a placeholder for N since we need 
% to change N for PotentialPenning(...) in this algorithm 

S = floor((sqrt(9-12*(1-Nf))-3)/6);           % Maximum number of complete subshells for Nf
[trash, x, y, N] = generate2DHexLattice(S,1); % Generate base lattice
Nleft = Nf - N;                               % Number of ions left to add
[trash x0 y0] = addHexShell(S+1);             % Generate next complete subshell (reffered to as current subshell)

for j=1:Nleft                % For each of the ions left to add
    V = zeros(1,length(x0)); % Create array to store potential energies
    N = N + 1;               % PotentialPenning will compute energy with one more ion 
    M = ones(1,N);           % Resize mass array

    for i=1:length(x0)       % For the points remaining the current subshell
        V(i) = PotentialPenning([x x0(i) y y0(i)]); % Find the potential energy associated with adding it
    end
    
	% Find the configuration with minimum potential energy 
    [Y,ind] = min(V); 

    % Permanently add this ion to the crystal
    x = [x x0(ind(1))];
    y = [y y0(ind(1))];
    
    % Remove that ion from the current subshell
    x0(ind(1))=[];
    y0(ind(1))=[];    
    
end

% Rescale entire configuration
u0 = [x y]*alpha;

end

