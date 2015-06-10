clf
setTrapParameters(0,0,0)
global G
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_16_NormalModeExpansionTest_inPlane\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_27_NormalModeExpansionTest_inPlane\';
FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_8_EnergyTest\';
%FileLocation = 'D:\PenningSimulationData\2014_3_28_SmallCrystalModes\';



%uavg = dlmread([FileLocation 'uavg.dat']);
Xpsd = dlmread([FileLocation 'Xpsd.dat']);
Ypsd = dlmread([FileLocation 'Ypsd.dat']);
params = dlmread([FileLocation 'params.dat']);
% params(1) = 127;
% params(2) = 44;
% params(3) = 0.0036;
% params(4) = 0;
% params(5) = 5e5;
% params(6) = 100;
% params(7) = 5e-9;


setTrapParameters(params(2),-params(3)/G,params(1));
%u0 = generateLattice(params(1),1);
%u = findEquilibrium(u0);
filename =[FileLocation int2str(params(5)-1) '.dat']; % use last configuration to find eigenvectors
%filename =[FileLocation int2str(66000) '.dat']; % use first configuration to find eigenvectors
M = dlmread(filename);
[u z] = convertPythonDataToMatlab(M);
u = findEquilibrium(u);
[E,D,st] = normalModes(u,0);

%uavg  = convertPythonDataToMatlab(uavg); % use average position for equillibrium positions
%[E,D,st] = normalModes(uavg,0);
global wz

freq = (0:0.5/(params(6)*params(7))/length(Xpsd):0.5/(params(6)*params(7)));
freq = 1.0e-6*freq(1:end-1);
figure(1)
semilogy(freq,sqrt(Xpsd.^2+Ypsd.^2),'k')
hold on
%semilogy(freq,Ypsd,'b')
for i=1:2*params(1)
    plot([1e-6*wz/(2*pi)*D(i) 1e-6*wz/(2*pi)*D(i)],[1e-15,1e3],'g')
end
%axis([0.6 .8 1e-11 1e-2])
xlabel('Frequency MHz','FontSize',24)
ylabel('Planar PSD','FontSize',24)
%title(['Magnetron Modes'],'FontSize',24)
title(['Cyclotron Modes'],'FontSize',24)
%title([' Axial Freq, N = ' num2str(params(1)) ', w = ' num2str(params(1)) ' kHz, delta = ' num2str(params(3))],'FontSize',24)
set(gca,'FontSize',24)

