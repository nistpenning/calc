% neighs = findNearestNeighbors(u)
% 
% Finds the nearest neighbors for each ion
%
% INPUT
% u:      Lattice
% 
% OUTPUT
% neighs: Cell list of neighboring indices of the ith ion
% 
% For example, neighs{1} gives the indices of u of neighboring particles of
% particle 1. That is, to access the positions of the neighboring
% particles, use x_neigh = u(neighs{i}) and y_neigh = u(neighs{i}+N) 

function neighs = findNearestNeighbors(u)

N=length(u)/2;                             %number of ions
neighs = cell(1,N);                        %Cell array/list of neighbors
tri = delaunay(u(1:end/2),u(end/2+1:end)); %find delaunay triangulation
% relevant MATLAB commands:
% 	triplot(tri,u(1:end/2),u(end/2+1:end))
% 	voronoi(x,y)

% tri is a mtri by 3 matrix, for mtri triangles with each row containing
% the indices of u that create that triangle

while ~isempty(tri)
    % for each line, store neighbors for each column
    neighs{tri(1,1)} = [neighs{tri(1,1)} tri(1,2:end)];
    neighs{tri(1,2)} = [neighs{tri(1,2)} tri(1,[1,3])];
    neighs{tri(1,3)} = [neighs{tri(1,3)} tri(1,1:2)];   
    tri(1,:)=[]; %delete current line
end

%because triangulation repeats indices, find unique neighboring indices
for i=1:N
    neighs{i} = unique( neighs{i} ); 
end


end