% [avg, R, avgdist]  = averageDistance(u)
%
% Finds mean distance of each ion from its nearest neighbors, then takes the
% average over all ions to find the mean separation, avgdist
%
% INPUT
% u:       Equilibrium lattice
%
% OUTPUT
% avg:     mean distance to neighbors for each ion
% R:       distance from origin for each ion
% avgdist: mean of avg

function [avg, R, avgdist] = averageDistance(u)

neighs = findNearestNeighbors(u);    % Find nearest neighbors
N = length(u)/2;                     % Number of ions
avg = zeros(1,N);                    % Stores mean distance to neighbors  

r = pairDistance(u);                 % Find pairwise distances
R = sqrt(u(1:N).^2 + u(N+1:end).^2); % Find distance from origin

for i=1:N
   avg(i)=mean(r(i,neighs{i}));      % Find mean distance from neighbors
end

[R,sorted] = sort(R);   			 % Sort R ascending
avg = avg(sorted);                   % Match avg with R
avgdist = mean(avg);    			 % Mean of avgdist (mean inter lattice spacing of equilbrium lattice)

end
