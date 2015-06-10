% z = getZforAxialMode(u,mode)
%
% Get z positions for axial mode <mode> at time zero, phase zero
% for dimensionless equilibrium positions u 
%
% mode can be vector of modes - will get superposition of modes for z
%
% Must have trap parameters set by setTrapParameters(...)

function z = getZforAxialMode(u,mode)

E = normalModes(u,1);
N = length(u)/2;

if (mode == 'all')
    mode=1:N;
end

zmax = 1e-8/length(mode);     % largest z value to have in crystal

z = zeros(N,1);

for m = mode
    %z = z + E(:,m)/max(abs(E(:,m)))*zmax; % scale 
    z = z + E(:,m)*zmax; % scale 
end
        
end



