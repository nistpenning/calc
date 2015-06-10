% Save configuration in data file in format that dominics code can read it
%
% Call setTrapParameters()

function pythonReadableU(u,z, saveas)

global N l0 q m

u = u*l0;

if z == 0
    z = zeros(1,N);
end

M = zeros(8,N);
M(1,:) = u(1:N);
M(2,:) = u(N+1:end);
M(3,:) = z; % should already be scaled
% same for 4-6 no momentum

M(7,:) = q*ones(1,N);
M(8,:) = m*ones(1,N);

dlmwrite(saveas,M,' ');

end

