setTrapParameters(0,0,0)
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_31_NonEquilibriumTest\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_10_AxialTemp_LaserCooling\';    

%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_23_NormalModeExpansionTest_BigData_smallZ-8\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_24_PlaneTransition\';
FileLocation = 'D:\PenningSimulationData\2014_3_28_SmallCrystalModes\';

psd = dlmread([FileLocation 'psd.dat']);
params = dlmread([FileLocation 'params.dat']);
% params(1) = 127;
% params(2) = 44;
% params(3) = 0.0036;
% params(4) = 0;
% params(5) = 5e5;
% params(6) = 100;
% params(7) = 5e-9;

global G

setTrapParameters(params(2),-params(3)/G,params(1));
%u0 = generateLattice(params(1),1);
%u = findEquilibrium(u0);
filename =[FileLocation int2str(params(5)-1) '.dat']; % use last configuration to find eigenvectors
M = dlmread(filename);
[u z] = convertPythonDataToMatlab(M);
[E,D,st] = normalModes(u,1);
global wz

freq = (0:0.5/(params(6)*params(7))/length(psd):0.5/(params(6)*params(7)));
freq = 1.0e-6*freq(1:end-1);
semilogy(freq,psd,'k')
hold on
for i=1:params(1)
    plot([1e-6*wz/(2*pi)*D(i) 1e-6*wz/(2*pi)*D(i)],[1e-16,1],'g')
end
%axis([0.6 .8 1e-11 1e-2])
axis([0.1 1 1e-16 1e-2])
xlabel('Frequency MHz','FontSize',24)
ylabel('PSD','FontSize',24)
title(['Compare Axial Freq, N = ' num2str(params(1)) ', w = ' num2str(params(2)) ' kHz, \delta = ' num2str(params(3))],'FontSize',24)
set(gca,'FontSize',24)

