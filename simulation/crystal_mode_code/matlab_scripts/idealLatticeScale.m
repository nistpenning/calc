% [u0 rescale] = idealLatticeScale(u0,min_scale,res)
%
% Ranges over several scaling constants, calculating the potential energy
% of that lattice (usually the seed lattice) with the trap parameters
% set by setTrapParameters. This helps the solver at calculating 
% the equilibrium structure, but it may not be entirely useful
% since the solver seems to be pretty robust.
%
% If a lattice constant was previously set, this function
% finds the scaling factor with respect to that lattice constant
%
% INPUT
% u0:        Seed lattice
% min_scale: Your guess at the possible minimum scaling constant
% res:       Resolution in searching for rescaling factor
%
% OUTPUT
% u0:        Rescaled lattice
% rescale:   Rescaling factor
%
% REMEMBER TO SET TRAP PARAMETERS (setTrapParameters(...))

function [u0 rescale] = idealLatticeScale(u0,min_scale,res)

for scale = min_scale:res:1000*min_scale % Scaling range
    next = PotentialPenning(u0*scale);   % Find potential energy of scaled lattice
    if scale~=min_scale                  % If not first iteration
        if next >= last                  % When next value is greater than before, minimum found
										 % since V is ~quadratic
            rescale = scale-res;         % Find lattice scaling, (last scaling)
            break; % FOUND
        end
    end
    last = next;   % Update last potential energy
end

u0 = u0*rescale;   % Rescale
end




