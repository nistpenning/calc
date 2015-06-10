% unew = rotate(u,theta)
%
% Rotates x, y coordinates by theta.
% Used to keep equilibrium positions on fixed orientation

function unew = rotate(u,theta)

x = u(1:end/2);
y = u(end/2+1:end);

xnew = x.*cos(theta) - y.*sin(theta);
ynew = x.*sin(theta) + y.*cos(theta);

unew = [xnew ynew];

end
