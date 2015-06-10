% Convert python data file and to format for matlab scripts
%
% Read parameters first and call setTrapParameters()

function [u,z,v] = convertPythonDataToMatlab(M)

global l0 v0

dim = size(M);
z = [];
v = [];

% Convert in plane positions
u = M(1,:);     % First row is x values
u = [u M(2,:)]; % Second row is y values, append to u
u = u/l0;       % Divide by characteristic length
if dim(1)>2
    z = M(3,:);    
    z = z/l0;
end

% if dim(1) > 3
% % Convert in plane velocities
% v = M(4,:);     % First row is vx values
% v = [v M(5,:)]; % Second row is vy values
% v = v/v0;       % Divide by characteristic velocity
% end

end

