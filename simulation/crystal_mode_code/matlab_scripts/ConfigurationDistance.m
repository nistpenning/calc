% ConfigurationDistance(u1,u2)
%
% Finds euclidean distance between two configurations
%
% Assumes u1 and u2 have same labeling of particles

function dist = ConfigurationDistance(u1,u2)

    dist = sqrt(sum((u2-u1).^2));

end
