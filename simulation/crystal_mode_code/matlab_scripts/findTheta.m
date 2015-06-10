% [theta u] = findTheta(u)
%
% Rotates the cloud to find the theta that gives the minimum potential
% energy. Used to fix unaligned clouds.
%
% INPUT
% u:        Input lattice
%
% OUTPUT
% theta:    rotate u0 by theta (clockwise) to give minimum potential
% u:        rotated crystal coordinates
%
% REMEMBER TO SET TRAP PARAMETERS (setTrapParameters(...))

function [theta u] = findTheta(u0)
                 
num = 200;                  % Number of points to scan
d_theta = pi/num;           % Resolution of rotation
potentials = zeros(num+1);    % Store energies
potentials(1) = PotentialPenning(u0); % Potential energy of unaligned crystal

for theta = d_theta:d_theta:pi                 % Scan through 0 to pi should because trap is symmetric along x and y
    
    u = rotate(u0,theta);              % Rotate original crystal by theta
    potentials(uint8(theta/d_theta)+1) = PotentialPenning(u); % save potential energy of this config

end

[trash ind] = min(potentials);   % find min potential
theta = (ind(1)-1)*d_theta;      % find theta for that potential
u = rotate(u0,theta);            % rotate to that config

end




