% Plot extracted amplitudes for each mode over the time bin
% to see how amplitude for each mode changes over time

setTrapParameters(0,0,0);
global m wz

hbar = 1.05457173e-34;
kB = 1.38e-23;
binsize = 1000;  

%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_5_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_9_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_5_NonEquilibriumTest_Mode45_zvels\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_6_NormalModeExpansionTest_zvels\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_10_AxialTemp_LaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_5_EnergyTest\';
%FileLocation = 'D:\PenningSimulationData\2014_3_27_AxialTemp_LaserCooling\';    
%FileLocation = 'D:\PenningSimulationData\2014_6_22_PlanarTemp\'; 
FileLocation = 'D:\PenningSimulationData\2014_7_15_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarTemp\'; 

params = dlmread([FileLocation 'params.dat']);
thetas = dlmread([FileLocation 'thetas.dat']);
%load([FileLocation 'axialModeDecomposition_eq.mat'],'norm_coords')
load([FileLocation 'axialModeDecomposition.mat'],'norm_coords')
global G N q V0
setTrapParameters(params(2),-params(3)/G,params(1));
% % If haven't calculated amps and freqs, do so now
% if ~exist([FileLocation 'amps' num2str(params(5)/binsize) '.dat'],'file')
%     disp('Calculating Amplitudes according to binsize')
%     axialModeAmplitude(FileLocation,1:params(1),binsize);
% end

%amps = dlmread([FileLocation 'amps' num2str(params(5)/binsize) '.dat']);
% freqs = dlmread([FileLocation 'freqs' num2str(params(5)/binsize) '.dat']);
% amps = dlmread([FileLocation 'maxamps' num2str(params(5)/binsize) '.dat']);
% 
% filename = [FileLocation int2str(params(5)-1) '.dat']; 
% M = dlmread(filename);
% u = convertPythonDataToMatlab(M);
% u = rotate(u,-thetas(params(5)-1)); 
% [E,D,st] = normalModes(u,1);

E = Ea;
D = Da;
modeTemp = [];
for j = 1:params(1)
    modeTemp(j) = 0.5*m*(wz*D(j))^2*mean(norm_coords(length(norm_coords)/10:end,j)'.^2)/(kB/2);
end

% 
% ionEnergy = zeros(1,params(5)/binsize); % Store energies calculated from ion velocities and potentials
% modeEnergy =  zeros(params(5)/binsize,127);
% totalModeEnergy = zeros(1,params(5)/binsize); % Store energies calculated from amplitude of harmonic oscillation from modes
% 
% for i = 1:params(5)/binsize
% %for i = 1:500000
% 
%     filename = [FileLocation int2str((i-1)*params(5)/binsize) '.dat']; 
%     %filename = [FileLocation int2str(i) '.dat']; 
% 
%     M = dlmread(filename);
%     z = M(3,:);
%     %vz = 0;
%     vz = M(4,:);
%     %energy(i,:) = ionAxialEnergy(z,vz);
%     
%     %ionEnergy(i) = sum(energy(i,:));
% %    ionEnergy(i)= 0.5*m*sum(vz.^2);
%  ionEnergy(i)= q*V0*sum(z.^2);
%     %scatter(i,ionEnergy/kB)
%     %  hold off
% % %     %plot(energy)
%     %  semilogy(energy)
%     %  hold on
%     
%     for j = 1:127
%         %const = 2*sqrt(hbar/(2*m*wz*freqs(i,j)));
%         const = 2*sqrt(hbar/(2*m*wz*D(j)));
%         %semilogy([j j],[1e-12 abs(amps(i,j))])
%         %plot([j j],[0 abs(amps(i,j))])
%         
%         %plot([wz*freqs(i,j)/2/pi wz*freqs(i,j)/2/pi],[0 abs(amps(i,j))])
%         %semilogy([wz*freqs(i,j)/2/pi wz*freqs(i,j)/2/pi],[1-12 abs(amps(i,j))])
%         %plot([j j],[0 abs(hbar*wz*freqs(j)*abs(amps(i,j))/const)])
%         %semilogy([j j],[1e-34 abs(hbar*wz*freqs(j)*abs(amps(i,j))/const)])
%         %semilogy([j j],[1e-34 abs(0.5*m*(wz*freqs(j))^2*abs(amps(i,j)).^2)])
%        % hold on
%         %E(i) = E(i) + hbar*wz*freqs(j)*abs(amps(i,j))/const;
%         %modeEnergy(i,j) = hbar*wz*freqs(i,j)*(abs(amps(i,j))/const).^2;
%         modeEnergy(i,j) = hbar*wz*D(j)*(abs(norm_coords(i,j)/const)).^2;
%         %modeEnergyClass(i,j) = 0.5*m*(wz*freqs(i,j))^2*abs(amps(i,j)).^2;
%         totalModeEnergy(i) = totalModeEnergy(i) + modeEnergy(i,j);
%         
%     end
%   % title([num2str(E(i)) ' ' num2str(ionEnergy(i))])
%     %axis([1 127 1e-12 1e-6])
%     %axis([1 127 0 5e-6])
%     %axis([1 127 1e-35 1e-22])
%     %axis([0.6e6 0.8e6  0 4e-6])
%     %pause(.000001)
%     %hold off
% end


