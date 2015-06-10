
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_16_AxialTemp_LaserCooling\';    
FileLocation = 'D:\PenningSimulationData\2014_3_27_AxialTemp_LaserCooling\';    
thetas = dlmread([FileLocation 'thetas.dat']);
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(0,0,0);
global G N
setTrapParameters(params(2),-params(3)/G,params(1));

norm_coords = zeros(params(5),params(1)); 
norm_vels = zeros(params(5),params(1)); 
% Use last configuration to find eigenvectors
filename = [FileLocation int2str(params(5)-1) '.dat']; 
M = dlmread(filename);
u = convertPythonDataToMatlab(M);
u = rotate(u,-thetas(params(5)-1)); 
u = findEquilibrium(u);                % original equilibrium lattice
[E,D,st] = normalModes(u,1);

for k = 1:params(5)
 
filename =[FileLocation int2str(k-1) '.dat'];
M = dlmread(filename);
z = M(3,:);    
v = M(6,:);    
norm_coords(k,:) = modeBasis(z,E);
norm_vels(k,:) = modeBasis(v,E);

end

save([FileLocation 'axialModeDecomposition_eq.mat'],'norm_coords')
save([FileLocation 'axialVelModeDecomposition_eq.mat'],'norm_vels')