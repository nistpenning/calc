% Compute equilibrium positions a lot, see the "landscape" of potential
% energies


%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_10_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_16_AxialTemp_LaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_8_EnergyTest\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_27_NormalModeExpansionTest_inPlane\';
%FileLocation = 'D:\PenningSimulationData\2014_3_25_AxialTemp_LaserCooling\';    
FileLocation = 'D:\PenningSimulationData\2014_3_28_SmallCrystalModes\';

thetas = dlmread([FileLocation 'thetas.dat']);
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(0,0,0);
global m wz l0 q V0 G 
setTrapParameters(params(2),-params(3)/G,params(1));

binsize = 1000;
pot = zeros(1,params(5)/binsize);
pot_ = zeros(1,params(5)/binsize);
eqdists = zeros(1,params(5)/binsize);

filename =[FileLocation int2str(0) '.dat'];
M = dlmread(filename);
u = convertPythonDataToMatlab(M);
u = rotate(u,-thetas(1));            % rotate to that config
u0 = findEquilibrium(u);                % original equilibrium lattice

for i = 0:binsize:params(5)-1

    filename =[FileLocation int2str(i) '.dat'];
    M = dlmread(filename);
    u = convertPythonDataToMatlab(M);
    u = rotate(u,-thetas(i+1));            % rotate to that config
    pot_((i/1000)+1) = PotentialPenning(u);
    [u pot((i/1000)+1)] = findEquilibrium(u);

    eqdists((i/1000)+1) = ConfigurationDistance(u0,u);


end

save([FileLocation 'LocalMinimums.mat'],'eqdists','pot')
plot(pot)



