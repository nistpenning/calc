% [r,xsquare,ysquare] = pairDistance(u)
%
% Calculates pairwise distances between ions
%
% INPUT
% u:       Lattice
% 
% OUTPUT
% r:       N x N matrix of pairwise distance.   r(i,j) = r(j,i) distance from ion i to ion j
% xsquare: N x N matrix of pairwise x distance. xsquare(i,j) = -xsquare(j,i)
% ysquare: N x N matrix of pairwise y distance. ysquare(i,j) = -ysquare(j,i)

function [r,xsquare,ysquare,x,y] = pairDistance(u)

N = length(u)/2;    % Number of ions
x = u(1:N);         % x positions
y = u(N+1:end);     % y positions

xsquare = zeros(N); % N x N pairwise x distance matrix
ysquare = zeros(N); % N x N pairwise y distance matrix

for i=1:N
    xsquare(i,:) = x(i) - x; % pairwise x distance between particles 
    ysquare(i,:) = y(i) - y; % pairwise y distance between particles 
end

r=sqrt(xsquare.^2 + ysquare.^2); %square matrix of radial distances

end