% [uadd x y] = addHexShell(s)
%
% Finds complete hex shell s
%
% x, y are lattice positions in each direction
% uadd is a single array of these positions with format
%      [x1 y1 x2 y2 ... xN yN]

function [uadd x y] = addHexShell(s)

a=0; b=0;  % coefficients of unit vectors
		   % a corresponds to unit vector (0,1)
		   % b corresponds to unit vector (sqrt(3)/2 , 1/2)
		   % (a, b) = (0, 0) is center point

% Simple pattern to define each shell
a = [s:-1:-s -s*ones(1,s-1) -s:s s*ones(1,s-1)];
b = [0:s s*ones(1,s-1) s:-1:-s -s*ones(1,s-1) -s:-1];

% Convert coefficients into position vectors using basis
x = (sqrt(3)/2).*b; 
y =  0.5.*b + a;   

% Convert to format to add next shell
uadd = zeros(1,2*length(a)); 
uadd(1:2:end) = x;  % Format: [x1 y1 x2 y2 ... xN yN]
uadd(2:2:end) = y;  % YES, this is the format for adding a complete hex shell, but not the final format

end