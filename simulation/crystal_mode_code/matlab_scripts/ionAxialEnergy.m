% energy = ionAxialEnergy(z,vz)
%  
% Calculate axial energy for each ion from trap potentials
%
% energy is a vector of energy for each ion
%
% Call setTrapParameters(...) first!

function energy = ionAxialEnergy(z,vz)

global m q V0

energy = q*V0*z.^2 + 0.5*m*vz.^2;
    
end