% [u V] = findEquilibrium(u0)
%
% Finds equilibrium lattice given seed lattice u0
% Calls "PotentialPenning" which must have had trap parameters 
% set by setTrapParameters(...)
%
% INPUT
% u0:        Seed lattice
%
% OUTPUT
% u:         Equilibrium Lattice
% V:         Potential energy of equilibrium configuration
% 
% You can play with solver options inside this function

function [u V] = findEquilibrium(u0)

N = length(u0)/2;          % Number of ions
%funtolerance = 1e-20;       % Tolerance of solver
funtolerance = 1e-19;       % Tolerance of solver

% Solver options
options = optimset('GradObj','on','Hessian','on',...
    'MaxFunEvals',100000*2*N,'MaxIter',10000,...
    'TolFun',funtolerance,'Tolx',funtolerance,'Display','off');

% MATLAB Solver "fminunc" finds local minimum
[u, V,exitflag,output,grad,hessian] = fminunc(@PotentialPenning,u0,options);
%exitflag
%output
norm(grad)


end