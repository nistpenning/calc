% Master Script
%% Load Data
cleanupObj = onCleanup(@CLEAN);

%FileLocation = 'D:\PenningSimulationData\2014_4_17_SmallCrystalModes\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_20_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_19_PlanarExcitation_3\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_22_PlanarTemp\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_15_PlanarTemp\';  
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarTemp\';  
FileLocation = 'D:\PenningSimulationData\2014_10_06_PlanarTemp\';  
setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);

%params(5)=100000;
%params(5)=200000;
v = zeros(params(5),2*params(1));

setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
disp('Loading Simulation Data...')
[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);

%[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,6);
kB = 1.38e-23;
disp('Finished Loading')
usunrotated = us;
% Move positions to rotating frame
V_0 = PotentialPenning( rotate(us(1,:),-thetas(1)));
for i = 0:params(5)-1
%for i = 0:200000
    %vlab(i+1,:) = [vxs(i+1,:) vys(i+1,:)]; % lab velocities
    %v(i+1,:) = spin_down(us(i+1,:),[vxs(i+1,:) vys(i+1,:)],params(2)); % spin down velocities
    us(i+1,:) = rotate(us(i+1,:),-thetas(i+1));           
    
    %PPE(i+1) = 0.5*m*wz^2*l0^2*(PotentialPenning( us(i+1,:)) - V_0)/kB;
    %R = sqrt(us(i+1,1:N).^2 + us(i+1,N+1:end).^2); % Find distance from origin
    %avgR(i+1) = mean(R);
    %[junk1, junk2, avgdist(i+1)] = averageDistance(us(i+1,:));
end
%PKE = 0.5*m*sum(v.^2,2)/kB;
%% Calculate Normal Modes for a particular Configuration
%u0 = findEquilibrium(us(1,:));
% Use first configuration, equilibrium positions
[Ea,Da] = normalModes(u0,1); % Axial Modes
[Ep,Dp] = normalModes(u0,0); % Planar Modes

%% View Simulation
disp('Viewing Data')
    setTrapParameters(54.5,-400,19);
    global N m wr wc q ke B wz V0 w G Vw ww M l0 v0
    u0 = generateLattice(N,1);
    u0 = rotate(u0,pi/2);     % get same major axis as dominics code
    u0 = findEquilibrium(u0);
    
figure
% vidpath = 'Heatup19';
%  vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
%  vwr.FrameRate = 30;
%  open(vwr);
%figure
%for i = 0
%for i = 0:1:params(5)-1
for i = 1:100:params(5)
timechunk = 200;
%for i = [1000:(1000+timechunk) 10000:(10000+timechunk) 100000:(100000+timechunk) 150000:(150000+timechunk) 160000:(160000+timechunk) 170000:(170000+timechunk) 180000:(180000+timechunk)]
%for i = 0:5:1000
    u = us(i+1,:);
    %u = rotate(us(i+1,:),-thetas(i+1));        
    scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
    hold on
    %scatter(u0(1:end/2),u0(end/2+1:end),'bo', 'filled');
    %quiver(u(1:N),u(N+1:2*N),v(i+1,1:N),v(i+1,1+N:end),0,'k','LineWidth',1)
    hold off
    %axis([-1.2*max(u),1.2*max(u),-1.2*max(u),1.2*max(u)])
    axis([-9,9,-9,9])
    %axis(15*[-1,1,-1,1])
    %axis([-2,2,-2,2])
    %axis([-4,4,-4,4])
    %axis(6*[-1,1,-1,1])
    title(['Rotating Frame, N = ' num2str(params(1)) ', \delta = ' num2str(params(3)) ', w = ' num2str(params(2)) ' kHz, i = ' num2str(i)])
    xlabel('Scaled position')
    %writeVideo(vwr,getframe(gcf));
    pause(.000001)
end
%close(vwr)

%%  Excite Tilt Mode for varying Vwall for Python Simulation 5/14/14

N = 7;
f = 64;
zmax = 5e-6;
date = '2014_5_14_TiltCouplingVwall/';
mode1 = 6;
for Vwall = [-20 -40 -60 -80 -100 -120 -140 -160]
    setTrapParameters(64,Vwall,N);
    global N m wr wc q ke B wz V0 w G Vw ww M l0 v0
    u0 = generateLattice(N,1);
    u0 = rotate(u0,pi/2);     % get same major axis as dominics code
    u0 = findEquilibrium(u0);
    [Ea,Da] = normalModes(u0,1); % Axial Modes
    filename = ['D:\PenningSimulationData\' date '\' num2str(N) '_' num2str(f) '_' num2str(Vwall) '_Amode' num2str(mode1) '.dat'];
    z = Ea(:,mode1)*zmax; % scale 
    pythonReadableU(u0,z,filename)
end

%%  Excite Tilt Mode for varying initial excitation for Python Simulation 5/23/14

N = 7;
f = 64;

date = '2014_6_10_TiltCouplingExcitation/';
mode1 = 6;
Vwall = -80;
%for zmax = [5e-6 1e-6 5e-7; 
ind = 1;
%for zmax = 5*logspace(-7,-6,10)
for zmax = linspace(2e-6,3e-6,10)
    setTrapParameters(64,Vwall,N);
    global N m wr wc q ke B wz V0 w G Vw ww M l0 v0
    u0 = generateLattice(N,1);
    u0 = rotate(u0,pi/2);     % get same major axis as dominics code
    u0 = findEquilibrium(u0);
    [Ea,Da] = normalModes(u0,1); % Axial Modes
    filename = ['D:\PenningSimulationData\' date '\' num2str(N) '_' num2str(f) '_' num2str(Vwall) '_zmax' num2str(ind) '_Amode' num2str(mode1) '.dat'];

    z = Ea(:,mode1)*zmax; % scale 
    pythonReadableU(u0,z,filename)
    ind = ind +1;
end

%% Tilt Coupling 7 ions with varying Vwall or initial excitation (2014-5-15)

%FileLocation = 'D:\PenningSimulationData\2014_5_23_TiltCouplingExcitation\';
FileLocation = 'D:\PenningSimulationData\2014_6_10_TiltCouplingExcitation\';

setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
thetas = dlmread([FileLocation 'thetas.dat']);

mode1 = 6;
Vwall=-80;
%for Vwall = [-20, -40, -60, -80, -100, -120, -140, -160]
for ind = 1:10

    setTrapParameters(params(2),Vwall,params(1));
    %[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,Vwall); %This function does not distinguish between text Data<>.dat
    [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,ind); %This function does not distinguish between text Data<>.dat
    u0 = findEquilibrium(us(1,:));
    [Ea,Da] = normalModes(u0,1); % Axial Modes
    
    %Project Axial Motion in to Normal Modes
    norm_coords = zeros(params(5),params(1)); 
    norm_vels = zeros(params(5),params(1)); 
    
    for k = 1:params(5)

        norm_coords(k,:) = modeBasis(zs(k,:),Ea);
        norm_vels(k,:) = modeBasis(vzs(k,:),Ea);

    end

    save([FileLocation 'axialModeDecomposition_amode6_Scale' num2str(ind) '.mat'],'norm_coords')
    save([FileLocation 'axialVelModeDecomposition_amode6_Scale' num2str(ind) '.mat'],'norm_vels')
    
%     save([FileLocation 'axialModeDecomposition_amode6_Vwall ' num2str(Vwall) '.mat'],'norm_coords')
%     save([FileLocation 'axialVelModeDecomposition_amode6_Vwall ' num2str(Vwall) '.mat'],'norm_vels')
    
    
end
   


%%  Excite Axial Modes for Python Simulation

N = 7;
f = 64;
Vwall = -80;
%setTrapParameters(54.5,-80,N);
%setTrapParameters(44,-80,N);
setTrapParameters(64,Vwall,N);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0
u0 = generateLattice(N,1);
u0 = rotate(u0,pi/2);     % get same major axis as dominics code
u0 = findEquilibrium(u0);
[Ea,Da] = normalModes(u0,1); % Axial Modes

date = '2014_5_14_SmallCrystalModes/';
%date = '2014_5_12_LargeCrystalModes/';

zmax = 5e-6;
for mode1 = 1:N
    %filename = ['D:\PenningSimulationData\' date '\19_54.5_-80_Amode' num2str(mode) '.dat'];
    %filename = ['D:\PenningSimulationData\' date '\' num2str(N) '_44_-80_Amode' num2str(mode) '.dat'];
    filename = ['D:\PenningSimulationData\' date '\' num2str(N) '_' num2str(f) '_' num2str(Vwall) '_Amode' num2str(mode1) '.dat'];
    z = Ea(:,mode1)*zmax; % scale 
    pythonReadableU(u0,z,filename)
end
        
%% Load Axial Excitations and plot PSD (2014-4-17)

FileLocation = 'D:\PenningSimulationData\2014_5_11_SmallCrystalModes\';  
setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);

freq = (0:0.5/(params(6)*params(7))/(params(5)/2):0.5/(params(6)*params(7)));
freq = 1.0e-6*freq(1:end-1); % chop of last one, because of Matlab ranges...
cmp = colormap(hsv(params(1)));
for mode1 = 1:params(1)
    mode1
    [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode1);
    u0 = findEquilibrium(us(1,:));
    [Ea,Da] = normalModes(u0,1); % Axial Modes
    
    %Project Axial Motion in to Normal Modes
    norm_coords = zeros(params(5),params(1)); 
    norm_vels = zeros(params(5),params(1)); 
    
    for k = 1:params(5)

        norm_coords(k,:) = modeBasis(zs(k,:),Ea);
        norm_vels(k,:) = modeBasis(vzs(k,:),Ea);

    end

    %save([FileLocation 'axialModeDecomposition_amode ' num2str(mode1) '.mat'],'norm_coords')
    %save([FileLocation 'axialVelModeDecomposition_amode ' num2str(mode1) '.mat'],'norm_vels')
    
    
    spectra = abs(fft(zs)).^2;
    Zpsd = sum(spectra, 2);
    Zpsd = Zpsd(1:(length(Zpsd)/2))+Zpsd(end:-1:length(Zpsd)/2+1);
    semilogy(freq,Zpsd,'Color',cmp(mode1,:))
    hold on
    pause(.01)
end
   
%     
for i=1:params(1)
    plot([1e-6*wz/(2*pi)*Da(i) 1e-6*wz/(2*pi)*Da(i)],[1e-20,1e-5],'k')
end


%% Find maximum amplitude over several periods for normal coordinates (also find energy)
setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
%FileLocation = 'D:\PenningSimulationData\2014_5_11_SmallCrystalModes\';
%FileLocation = 'D:\PenningSimulationData\2014_5_14_TiltCouplingVwall\';  

%FileLocation = 'D:\PenningSimulationData\2014_5_23_TiltCouplingExcitation\';
FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation\';  

%FileLocation = 'D:\PenningSimulationData\2014_6_10_TiltCouplingExcitation\';


params = dlmread([FileLocation 'params.dat']);
thetas = dlmread([FileLocation 'thetas.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));

%for mode = 1:params(1)
%for mode = 1:2*params(1)
for mode = 1
%mode1 = 6;
Das = [];
%for Vwall = [-20, -40, -60, -80, -100, -120, -140, -160]
%Vwall = -80;
%for scale = 1:10

%setTrapParameters(params(2),Vwall,params(1));

% Get initial crystal config
%filename = [FileLocation num2str(params(1)) '_' num2str(params(2)) '_' num2str(-params(3)/G) '_Amode' num2str(mode) '.dat'];
%filename = [FileLocation num2str(params(1)) '_' num2str(params(2)) '_' num2str(Vwall) '_Amode' num2str(mode1) '.dat'];
%filename = [FileLocation num2str(params(1)) '_' num2str(params(2)) '_' num2str(Vwall) '_zmax' num2str(scale) '_Amode' num2str(mode1) '.dat'];

% Mat = dlmread(filename);
% u = Mat(1,:);     % First row is x values
% u = [u Mat(2,:)]; % Second row is y values, append to u
% u0 = u/l0;        % Divide by characteristic length
% [Ea,Da] = normalModes(u0,1); % Axial Modes
% Das = [Das Da];
% tiltsep = wz*(Da(mode1) - Da(mode1-1))/2/pi;

% Load Normal Coordinates 
%load([FileLocation 'axialModeDecomposition_amode ' num2str(mode) '.mat'],'norm_coords')
%load([FileLocation 'axialModeDecomposition_amode6_Vwall ' num2str(Vwall) '.mat'],'norm_coords')
%load([FileLocation 'axialModeDecomposition_amode6_Scale' num2str(scale) '.mat'],'norm_coords')
%load([FileLocation 'axialModeDecomposition_amode6_Scale' num2str(scale) '.mat'],'norm_coords')
load([FileLocation 'axialMode_' num2str(mode) '.mat'],'norm_coords')

dims = size(norm_coords);

% scale = params(end);
% if abs((norm_coords(1,mode1)/scale)-1) > 1e-4;
%     error('Scaling does not match')
% end

ModeEnergies =  0.5*m*(norm_coords).^2.*repmat(((wz*Da).^2)',dims(1),1);
%E0_all = sum(ModeEnergies(1,:));
%E0 = ModeEnergies(1,mode); 
cutoff = 1e5; % 100 kHz
%ModeEnergyEnvelope = ModeEnergies/E0;
%ModeEnergyEnvelope = (2/E0)*filter(params(7)*params(6)*cutoff, [1 params(7)*params(6)*cutoff-1], ModeEnergies); % low pass filter, multiply by two to fix attenuation  
ModeEnergyEnvelope = 2*filter(params(7)*params(6)*cutoff, [1 params(7)*params(6)*cutoff-1], ModeEnergies); % low pass filter, multiply by two to fix attenuation  
TotalModeEnergy = sum(ModeEnergyEnvelope,2);

figure
cmp = colormap(jet(dims(2)));
set(gcf,'DefaultAxesColorOrder',cmp)
times = linspace(0,params(5)*params(7)*params(6),params(5));
semilogy(times, ModeEnergyEnvelope)
%semilogy(times, ModeEnergyEnvelope(:,[mode-1 mode]))
%axis([0,times(end),1e-10,1])

xlabel('Time (seconds)','FontSize',24)
ylabel('Energy Relative to Initial Excited Mode','FontSize',24)
%title(['Excited Axial Mode ' num2str(mode) ' of ' num2str(params(1)) ' by ' num2str(scale/1e-9) ' nm'], 'FontSize',14)
%title(['Tilt Excitation by ' num2str(scale/1e-9) ' nm at Tilt Frequency Seperation ' num2str(tiltsep) ' Hz'], 'FontSize',14)
test =  5*logspace(-7,-6,10);
%test = linspace(2e-6,3e-6,10);
%title(['Tilt Excitation by ' num2str(test(scale)/1e-9) ' nm'], 'FontSize',14)
clb = colorbar;
set(clb,'YTick',1:params(1))
set(gca,'FontSize',24)
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
%SaveLocation = 'C:\Users\ACKWinDesk\Google Drive\PenningSimulation\2014_05_14\';
%SaveLocation = 'C:\Users\ACKWinDesk\Google Drive\PenningSimulation\2014_06_10\';
%saveas(gcf,[SaveLocation 'DrivingTiltMode6_Scale_' num2str(round(test(scale)/1e-9))],'png')
%saveas(gcf,[SaveLocation 'DrivingEachMode' num2str(mode) 'of' num2str(params(1))],'png')
%saveas(gcf,[SaveLocation 'DrivingTiltMode6_Vwall' num2str(Vwall)],'png')
%pause(2)
%close(gcf)

end

% % Load crystal data for Planar Energy
% [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode);
% 
% PlanarKineticIons_rotframe = [];
% R = [];
% for i = 1:1000:params(5)
%     u = us(i,:);  % don't rotate, want origin lattice to spin down
%     %u = rotate(us(i,:),-thetas(i)); 
%     R = [R; sqrt(u(1:N).^2 + u(N+1:end).^2)]; % Find distance from origin
%     v = spin_down(u,[vxs(i,:) vys(i,:)],params(2));
%     vx = v(1:N);
%     vy = v(N+1:end);
%     PlanarKineticIons_rotframe = [PlanarKineticIons_rotframe; 0.5*m*(vx.^2+vy.^2)];
% end
% TotalPlanarKE = sum(PlanarKineticIons_rotframe,2)/E0;
% 
% 
% figure
% times2 = times(1:1000:end);
% semilogy(times2,TotalPlanarKE)
% hold on
% semilogy(times,TotalModeEnergy,'g')


% figure
% envelope(:,mode)=[];
% set(gcf,'DefaultAxesColorOrder',cmp)
% plot(times,envelope)
% set(gca,'FontSize',24)
% xlabel('Time (seconds)','FontSize',24)
% ylabel('Amplitude (scaled by initial scaling)','FontSize',24)
% title(['Excite Mode ' num2str(mode) ' each Axial mode by 10 nm'], 'FontSize',14)
% set(gca,'FontSize',24)

%E0 =  0.5*m*sum( ((wz*Da).^2).*(norm_coords(1,:)'.^2));
%envelope = reshape(max(reshape(norm_coords,blocksize,dims(1)*dims(2)/blocksize)),dims(1)/blocksize,dims(2))/scale; %find max amplitude over blocksize for each mode
%ModeEnergies =  0.5*m*(envelope*scale).^2.*repmat(((wz*Da).^2)',dims(1)/blocksize,1);
%TotalModeEnergy = sum(ModeEnergies,2)/E0;
%figure
%semilogy(TotalModeEnergy)

% blocksize = 1;
% if blocksize ~= 1
%     ModeEnergyEnvelope = reshape(max(reshape(ModeEnergies,blocksize,dims(1)*dims(2)/blocksize)),dims(1)/blocksize,dims(2))/E0; 
% else
%     ModeEnergyEnvelope = ModeEnergies/E0;
% end
% TotalModeEnergy = sum(ModeEnergyEnvelope,2);
%% Find maximum amplitude over several periods for planar normal coordinates (also find energy)
setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
%FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_29_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_13_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_9_24_PlanarExcitation\';  
FileLocation = 'D:\PenningSimulationData\2014_9_24_PlanarExcitation\';  
params = dlmread([FileLocation 'params.dat']);
hbar = 1.05457173e-34;
Tdoppler = hbar*2*pi*19e6/2/kB;
thetas = dlmread([FileLocation 'thetas.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));

%for mode = 1:2*params(1)
%for mode = 1:params(1)
%for mode = params(1)+1:2*params(1)

for mode = 4+N

Das = [];

%load([FileLocation 'planarMode_' num2str(mode) 'old.mat'],'norm_coords_planar','EnergyMode')
%load([FileLocation 'planarMode_' num2str(mode)'.mat'],'norm_coords_planar')
load([FileLocation 'planarMode_' num2str(mode) '.mat'],'norm_coords_planar','EnergyMode')
%dims = size(norm_coords);

% scale = params(end);
% if abs((norm_coords(1,mode1)/scale)-1) > 1e-4;
%     error('Scaling does not match')
% end

% EnergyMode = zeros(params(5),2*params(1));
% for k = 1:N
%     EnergyMode(:,k) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2);
%     %EnergyMode(:,k) = (abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2)*( (real((Ep(1:2*N,2*k-1))'*(Dp(2*k)*Dp(2*k))*(Ep(1:2*N,2*k-1))  -0.5*Ep(1:2*N,2*k-1)'*V*Ep(1:2*N,2*k-1))) + ...
%     %        (real((Ep(1:2*N,2*k))'*(Dp(2*k)*D(2*k))*(Ep(1:2*N,2*k))  -0.5*Ep(1:2*N,2*k)'*V*Ep(1:2*N,2*k))));
% end
EnergyMode(10000,mode)/Tdoppler

%E0 = ModeEnergies(1,mode); 
cutoff = 1e5; % 100 kHz
%cutoff = 1e9; % 100 kHz
%ModeEnergyEnvelope = ModeEnergies/E0;
%ModeEnergyEnvelope = (2/E0)*filter(params(7)*params(6)*cutoff, [1 params(7)*params(6)*cutoff-1], ModeEnergies); % low pass filter, multiply by two to fix attenuation  
ModeEnergyEnvelope = 2*filter(params(7)*params(6)*cutoff, [1 params(7)*params(6)*cutoff-1], EnergyMode/Tdoppler); % low pass filter, multiply by two to fix attenuation  
TotalModeEnergy = sum(ModeEnergyEnvelope,2);

figure
%cmp = colormap(jet(dims(2)));
%cmp = colormap(jet(2*N));
cmp = colormap(jet(N));
set(gcf,'DefaultAxesColorOrder',cmp)
times = linspace(0,params(5)*params(7)*params(6),params(5));
semilogy(times, ModeEnergyEnvelope)%(:,(2:end))
%loglog(times, ModeEnergyEnvelope)%(:,(2:end))
%semilogy(times, ModeEnergyEnvelope(:,[mode-1 mode]))
%axis([0,times(end),1e-10,1e5])
axis([0,times(end),1e-5,1e10])

xlabel('Time (seconds)','FontSize',24)
%ylabel('Energy Relative to Initial Excited Mode','FontSize',24)
ylabel('Energy (units of T_{Doppler})','FontSize',24)
%title(['Excited Planar Mode ' num2str(mode) ' of ' num2str(2*params(1))], 'FontSize',14)
title(['Excited Magnetron Mode ' num2str(mode) ' of ' num2str(params(1)) ', scale = 1e-0'], 'FontSize',14)

test =  5*logspace(-7,-6,10);
%test = linspace(2e-6,3e-6,10);
%title(['Tilt Excitation by ' num2str(test(scale)/1e-9) ' nm'], 'FontSize',14)
clb = colorbar;
%set(clb,'YTick',1:2*params(1))
set(clb,'YTick',1:params(1))
set(gca,'FontSize',24)
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

%SaveLocation = 'C:\Users\ACKWinDesk\Google Drive\PenningSimulation\2014_09_26\';
%saveas(gcf,[SaveLocation 'DrivingTiltMode6_Scale_' num2str(round(test(scale)/1e-9))],'png')
%saveas(gcf,[SaveLocation 'ExcitingEachMagnetronMode1e-0_' num2str(mode) 'of' num2str(params(1)) ],'png')
%saveas(gcf,[SaveLocation 'DrivingTiltMode6_Vwall' num2str(Vwall)],'png')
%pause(2)
%close(gcf)


end

%% Movie of Planar Modes

mg = abs(Ep);
ph = angle(Ep);

%vidpath = 'CyclotronModes19';
 %vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
 %vwr.FrameRate = 30;
 %open(vwr);

%for mode = N+1:2*N
for mode1 = 2
%u = us(1,:);
stretch = 0.5;
for wt = 0 : 2*pi/50 : 6*pi
    scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
    hold on
    scatter(u0(1:N)'+stretch*mg(1:N,mode1).*cos(ph(1:N,mode1)+wt),u0(N+1:2*N)'+stretch*mg(N+1:2*N,mode1).*cos(ph(N+1:2*N,mode1)+wt),'go', 'filled');

    hold off
    %quiver(u(1:N)',u(N+1:2*N)',stretch*mg(1:N,mode).*cos(ph(1:N,mode)+wt),...
    %    stretch*mg(N+1:2*N,mode).*cos(ph(N+1:2*N,mode)+wt),0,'k','LineWidth',1)
        axis([-1.2*max(u0),1.2*max(u0),-1.2*max(u0),1.2*max(u0)])
        title(['Cyclotron Mode: ' num2str(mode1)])
 %          writeVideo(vwr,getframe(gcf));
    pause(.001)
end
end
%close(vwr)


%% Radial Kinetic Energy vs Radius (new data loading)

PlanarKineticIons_rotframe = [];
R = [];

for i = 1:100:params(5)

    u = rotate(us(i,:),thetas(i));  % rotate back to spin down
    R = [R; sqrt(u(1:N).^2 + u(N+1:end).^2)]; % Find distance from origin
    v = spin_down(u,[vxs(i,:) vys(i,:)],params(2));
    vx = v(1:N);
    vy = v(N+1:end);
    PlanarKineticIons_rotframe = [PlanarKineticIons_rotframe; 0.5*m*(vx.^2+vy.^2)];
end

kB = 1.3806488e-23;

[Rsort ind] = sort(mean(R));
PlanarKinetic = mean(PlanarKineticIons_rotframe);
scatter(Rsort,PlanarKinetic/kB/1e-3)

%% Calculate PSDs
disp('Calculating PSDs')
freq = (0:0.5/(params(6)*params(7))/(params(5)/2):0.5/(params(6)*params(7)));
freq = 1.0e-6*freq(1:end-1); % chop of last one, because of Matlab ranges...

% Calculate PSD for Axial Motion
spectra = abs(fft(zs)).^2;
Zpsd = sum(spectra, 2);
Zpsd = Zpsd(1:(length(Zpsd)/2))+Zpsd(end:-1:length(Zpsd)/2+1);

% Calculate PSD for Planar Motion
motion = us - repmat(us(1,:),params(5),1); % subtract off equilibrium positions
spectra = abs(fft(motion)).^2;
Ppsd = sum(spectra, 2);
Ppsd = Ppsd(1:(length(Ppsd)/2))+Ppsd(end:-1:length(Ppsd)/2+1);
disp('Done Calculating PSDs')

%% Plot PSDs

% figure
% semilogy(freq,Zpsd)
% hold on
% for i=1:params(1)
%     plot([1e-6*wz/(2*pi)*Da(i) 1e-6*wz/(2*pi)*Da(i)],[1e-16,1],'g')
% end

figure
semilogy(freq,Ppsd)
hold on

for i=1:2*params(1)
    plot([1e-6*wz/(2*pi)*Dp(i) 1e-6*wz/(2*pi)*Dp(i)],[1e-15,1e15],'g')
end
xlabel('Frequency MHz','FontSize',24)
ylabel('Planar PSD','FontSize',24)
set(gca,'FontSize',24)
    

%% Load Data through files - for some reason this doesn't work very well
clear all

FileLocation = 'D:\PenningSimulationData\2014_3_28_SmallCrystalModes\';  
setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
disp('Loading Simulation Data')

us = zeros(params(5),2*params(1));
for i = 0:params(5)-1
    if ~mod(i,1000)
       disp(i)
    end
    filename =[FileLocation int2str(i) '.dat'];
    M = dlmread(filename);
    us(i+1,:) = convertPythonDataToMatlab(M);
    us(i+1,:) = rotate(us(i+1,:),-thetas(i+1));     
    us(i+1,:) = us(i+1,:) - us(1,:);
end

%% Radial Kinetic Energy vs Radius
clear all
setTrapParameters(0,0,0);
global G
FileLocation = 'D:\PenningSimulationData\2014_3_25_AxialTemp_LaserCooling\';   
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
N = params(1);
global m
%PlanarKineticIons_rotframe = zeros(params(5),N);
PlanarKineticIons_rotframe = [];
R = [];
for i = 1:1000:params(5)
    
    filename =[FileLocation int2str(i-1) '.dat'];
    M = dlmread(filename);
    u = convertPythonDataToMatlab(M); % dimensionless positions
    %u = rotate(u,-thetas(i)); 
    R = [R; sqrt(u(1:N).^2 + u(N+1:end).^2)]; % Find distance from origin

    vx = M(4,:);
    vy = M(5,:);
    v = spin_down(u,[vx vy],params(2));
    vx = v(1:N);
    vy = v(N+1:end);
    
    %PlanarKineticIons_rotframe(i,:) = 0.5*m*(vx.^2+vy.^2);
    PlanarKineticIons_rotframe = [PlanarKineticIons_rotframe; 0.5*m*(vx.^2+vy.^2)];
    
    
    if ~mod(i-1,1000)
       disp(i)
    end

end

[Rsort ind] = sort(mean(R));
PlanarKinetic = mean(PlanarKineticIons_rotframe);
scatter(Rsort,PlanarKinetic/kB/1e-3)

%% Calculate Equilibrium Positions and Save for Python Simulation

setTrapParameters(54.5,-80,19);
global N
u0 = generateLattice(N,1);
u0 = rotate(u0,pi/2);     % get same major axis as dominics code
u0 = findEquilibrium(u0);
filename = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\MatlabData\19_54.5_-80_eq.dat';
pythonReadableU(u0,0,filename)


%% Compare PolyEig to Matrix Pencil Linearization

N = 7;
%setTrapParameters(54.5,-80,N);
setTrapParameters(64,-80,N);
global N M 
u0 = generateLattice(N,1);
u0 = rotate(u0,pi/2);     % get same major axis as dominics code
u = findEquilibrium(u0);
%[Ep,Dp,st,A,T] = normalModes(u,0); % Planar Modes
[Ep,Dp,st,A,T] = normalModes_5_27_14(u,0); % Planar Modes
%pencil = [1i*T A; -eye(2*N) zeros(2*N)];
%A = A/2;
%pencil1 = [zeros(2*N) eye(2*N); A/2 -T]; %[E,D] = polyeig(-A/2,1i*T,-diag([M M]));
pencil1 = [-1i*T/2 -0.25*T*T+A/2; eye(2*N) -1i*T/2];
%pencil2 = [zeros(2*N) sqrt(A); -sqrt(A) T];
[E1, D1] = eig(pencil1);
%[E2, D2] = eig(pencil2);

D1 = diag(D1);
%return
D1 = real(D1);
%D1 = imag(D1);
gz = (D1>0);
E1 = E1(:,gz);
D1 = D1(gz);
[D1 ind] = sort(D1);
E1 = E1(:,ind); %sorted eigenvectors
Mass = diag([M M]);
disp('------------')
mode2=1;
mode1=1;
Energyp = imag(Ep(:,mode1))'*(Mass/(Dp(mode1)*Dp(mode2)))*imag(Ep(:,mode2)) - 0.5*real(Ep(1:2*N,mode1))'*A*real(Ep(1:2*N,mode2));
for mode1 = 1:2*N
    %mode2=mode1;
Energy1 = E1(1:2*N,mode1)'*(Mass/D1(mode1)^2 - A)*E1(1:2*N,mode1)/4;
%Energyp = imag(Ep(:,mode))'*(Mass/Dp(mode)^2)*imag(Ep(:,mode)) - 0.5*real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode));
Energyp = imag(Ep(:,mode1))'*(Mass/(Dp(mode1)*Dp(mode2)))*imag(Ep(:,mode2)) - 0.5*real(Ep(1:2*N,mode1))'*A*real(Ep(1:2*N,mode2));
%Energyp = 1;
%EnergyR = real(E1(2*N+1:end,mode))'*Mass*real(E1(2*N+1:end,mode)) - 0.5*real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode));
%EnergyI = 0.5*imag(E1(2*N+1:end,mode))'*Mass*imag(E1(2*N+1:end,mode)) - 0.5*imag(E1(1:2*N,mode))'*A*imag(E1(1:2*N,mode));

%Energy = E1(2*N+1:end,mode)'*Mass*E1(2*N+1:end,mode)/4;
%disp(['Energy = ' num2str(Energy)])
%disp(2*real(E1(2*N+1:end,mode))'*real(E1(2*N+1:end,mode))+real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode)))
%disp(abs(E1(2*N+1:end,mode))'*abs(E1(2*N+1:end,mode))+abs(E1(1:2*N,mode))'*A*abs(E1(1:2*N,mode)))

%disp((imag(E1(1:2*N,mode))'*imag(E1(1:2*N,mode))/D1(mode)^2 + real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode))))
%disp(['Overlap = ' num2str((imag(E1(1:2*N,mode))'*imag(E1(1:2*N,mode))/D1(mode)^2 - 0.5*real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode)))/(Energyp))])
%disp(['Overlap = ' num2str((real(E1(2*N+1:end,mode))'*real(E1(2*N+1:end,mode)) + real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode))))])

disp(['Overlap = ' num2str((imag(Ep(:,mode1))'*imag(Ep(:,mode2))./(Dp(mode1)*Dp(mode2)) - 0.5*real(Ep(:,mode1))'*A*real(Ep(:,mode2)))/Energyp)])

%real(E1(1:2*N,mode))'*A*real(E1(1:2*N,mode))
%real(E1(2*N+1:end,mode))'*real(E1(2*N+1:end,mode))
end



% 
% i = 1;
% j = 1;
% Mass = diag(M);
% real(E1(N+1:end,i))*Mass*real(E1(N+1:end,j)) + real(E1(1:N,i))*A*real(E1(1:N,j))


% D = diag(D);
% D = real(D);
% gz = (D>0);
% E = E(:,gz);
% D = D(gz);
% [D ind] = sort(D);
% E = E(:,ind); %sorted eigenvectors
% 
% ortho1 = zeros(2*N);
% ortho2 = zeros(2*N);
% for i=1:2*N
%     for j=1:2*N
%         ortho1(i,j) = dot(Ep(:,i),Ep(:,j));
%         ortho2(i,j) = dot(E(1:38,i),E(1:38,j));
%     end
% end

% real(E(1:2*N,mode))'*A*real(E(1:2*N,mode))
% real(E(2*N+1:end,mode))'*real(E(2*N+1:end,mode))
% 
% real(Ep(:,mode))'*A*real(Ep(:,mode))
% imag(Ep(:,mode))'*imag(Ep(:,mode))/Dp(mode)^2

% test = -eig(A);
% 
% l=1;
% k=2;
% %sum((Dp(l)*Dp(k) + test'.^2).*(conj(Ep(:,l)).*Ep(:,k)))
% ortho = 0;
% for i=1:2*19
%    ortho = ortho + (Dp(l)*Dp(k)+ test(i)^2)*(conj(Ep(i,l))*Ep(i,k));
% end

%% Movie of Axial Modes

figure;
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cmp = jet;
%set(gcf, 'Position',[.25*scrsz(3) 0 518   688]); % Maximize figure.
set(gcf,'color','k')

vidpath = 'Axial19';
vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
setTrapParameters(54.5,-400,19);
u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium
[Ea,Da] = normalModes(u0,1);

open(vwr);
%vwr.Quality = 100;
stretch = 50;
index = 1;   

for wt= 0 : 2*pi/50 : 2*pi
    ha = tight_subplot(4,5,[.01 .03],[.1 .01],[.01 .01]);
set(findobj(gcf, 'type','axes'), 'Visible','off')     
    iter = 1;
    %for i = [1 2 3 N-8:N]
    for i = 1:19
        z = Ea(:,i)'.*sin(wt);
        blu = max(abs(Ea(:,i)));
        
        %Color map HSV has 64 unique colors
        c = ceil(32*(z/blu +1));
        c(c<=1)=1;
        c(c>=64)=64;
        %cur = subplot(4,3,iter);
        %set(cur,'Position',[(0.05 + (mod(iter-1,3))*0.3), (0.05 + (2-floor(iter-1/3))*0.225), 0.28, 0.205])
        %axes('Position',[(0.001 + (mod(iter-1,3))*0.32), (0 + (5-floor((iter-1)/3))*0.2475), 0.315, 0.30])
        %axes('Position',[(0.001 + (mod(iter-1,3))*0.32), (0 + (5-floor((iter-1)/3))*0.2475), 0.315, 0.30])
        %subplot(4,5,i)
        axes(ha(i));
        scatter(u0(1:N),u0(N+1:2*N),(stretch/blu)*(blu+z),cmp(c,:),'filled');
       
        %hold on
        %scatter(0,0,'k+');
        %drawnow
        %hold off
        %axis tight
        %axis equal
        axis([-1.2*max(u0),1.2*max(u0),-1.2*max(u0),1.2*max(u0)])
        %out = 30;
        %axis([-out,out,-out,out])
        %axis equal
        %set(gca,'YTick',[])
        %set(gca,'XTick',[])
        %set(gca,'box','off')
        %axis off;

        %set(get(cur,'Parent'),'Position',[(0.05 + (mod(iter-1,3))*0.3), (0.05 + (2-floor(iter-1/3))*0.225), 0.28, 0.205])
        
        iter = iter+1;
        
    end
    pause(0.0001)
    
     set(findobj(gcf, 'type','axes'), 'Visible','off')     
    writeVideo(vwr,getframe(gcf));
    clf;
    index = index +1;
end
 %imwrite(im,map,'AxialTest.gif','DelayTime',0,'LoopCount',inf) %g443800
close(vwr);

%% Movie of Planar Modes

N = 7;  % Number of ions
f = 64;
Vwall = -80;
setTrapParameters(f,Vwall,N);% Set global variables
global M l0 q m
Mass = diag([M M]);

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium


[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction


vidpath = 'Magnetron7';
%vidpath = 'Cyclotron19';
vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
%vwr.FrameRate = 30;
open(vwr);
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w')
ha = tight_subplot(2,4,[.01 .03],[.1 .01],[.01 .01]);
%for mode = N+1:2*N
%scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
%hold on
for wt = 0 : 2*pi/50 : 2*pi
    for mode = 1:N
        %realE = 10*real(Ep(1:2*N,mode+N));
        %imagE = 10*imag(Ep(1:2*N,mode+N));
        realE = .5*real(Ep(1:2*N,mode));
        imagE = .5*imag(Ep(1:2*N,mode));
        %subplot(4,5,mode)
        axes(ha(mode));
        %axes('Position',[(0.001 + (mod(iter-1,3))*0.32), (0 + (5-floor((iter-1)/3))*0.2475), 0.315, 0.30])
        scatter(u0(1:N)'+ realE(1:N)*cos(wt) - imagE(1:N)*sin(wt),u0(N+1:2*N)'+ realE(N+1:2*N)*cos(wt) - imagE(N+1:2*N)*sin(wt),'bo', 'filled');
        %axis([-1.2*max(u0),1.2*max(u0),-1.2*max(u0),1.2*max(u0)])
        axis(1.5*[-max(u0),max(u0),-max(u0),max(u0)])
        %hold on
        %scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
    
    %title(['Magnetron Mode: ' num2str(mode1) ', Frequency = ' num2str(Dp(mode1)*wz/2/pi/1000) ' kHz, Rotation Frequency = 64 kHz, \delta = 0.0036' ])
    %title(['Cyclotron Mode: ' num2str(mode1-N) ', Frequency = ' num2str(Dp(mode1)*wz/2/pi/1000) ' kHz, Rotation Frequency = 64 kHz, \delta = 0.0036' ])
    %quiver(u(1:N)',u(N+1:2*N)',stretch*mg(1:N,mode).*cos(ph(1:N,mode)+wt),...
    %    stretch*mg(N+1:2*N,mode).*cos(ph(N+1:2*N,mode)+wt),0,'k','LineWidth',1)
        
      %  title(['Cyclotron Mode: ' num2str(mode1)])
    %writeVideo(vwr,getframe(gcf));
    %pause(.000001)
    end
    set(findobj(gcf, 'type','axes'), 'Visible','off')
    writeVideo(vwr,getframe(gcf));
    pause(0.0001)
end
close(vwr)
%pause(3)
%clf

%% Movie of Planar Modes

%N = 7;  % Number of ions
%f = 64;
%Vwall = -80;

N = 19;  % Number of ions
f = 54.5;
Vwall = -400;
%setTrapParameters(f,Vwall,N);% Set global variables
%N = 127;  % Number of ions
%setTrapParameters(44,-80,N);% Set global variables
global M l0 q m
Mass = diag([M M]);

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium

%[Eold,Dold,st,V,A] = normalModes_OLD(u0,0); % Planar Modes with correction
%[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes with correction
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction

%Ep = Em;
%mg = abs(E1);
%ph = angle(E1);

%vidpath = 'Magnetron19';
%vidpath = 'Cyclotron7';
%vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
%vwr.FrameRate = 30;
%open(vwr);
figure
%for mode = N+1:2*N
%scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
%hold on
%for mode1 = N+1:2*N
for mode1 = 1:2*N
%for mode1 = 28

scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
hold on
%for i = 1:N
    
realE = real(Ep(1:2*N,mode1));
imagE = imag(Ep(1:2*N,mode1));

% realEM = real(Eold(1:2*N,mode1));
% imagEM = imag(Eold(1:2*N,mode1));



%realEM = real(EM(1:2*N,mode1));
%imagEM = imag(EM(1:2*N,mode1));

%u = us(1,:);


%wt = 0 : 2*pi/100 : 2*pi;
%plot(u(i)+ realE(i).*cos(wt) - imagE(i).*sin(wt),u(i+N)+ realE(i+N).*cos(wt) - imagE(i+N).*sin(wt),'g-');
%plot(u(i)+ 5*realE(i).*cos(wt) - 5*imagE(i).*sin(wt),u(i+N)+ 5*realE(i+N).*cos(wt) - 5*imagE(i+N).*sin(wt),'g-');
%scatter(u(i)+ realE(i),u(i+N)+ realE(i+N),'bo','filled');
%quiver(u(i),u(i+N),5*realE(i), 5*realE(i+N),'b');
%axis([-1.2*max(u0),1.2*max(u0),-1.2*max(u0),1.2*max(u0)])
 
% stretch = 1;
for wt = 0 : 2*pi/50 : 100*pi
    figure(1)
    scatter(u0(1:end/2),u0(end/2+1:end),'ko', 'filled');
    hold on
    %scatter(u0(1:N)'+ 5*realE(1:N)*cos(wt) - 5*imagE(1:N)*sin(wt),u0(N+1:2*N)'+ 5*realE(N+1:2*N)*cos(wt) - 5*imagE(N+1:2*N)*sin(wt),'bo', 'filled');
    scatter(u0(1:N)'+ realE(1:N)*cos(wt) - imagE(1:N)*sin(wt),u0(N+1:2*N)'+ realE(N+1:2*N)*cos(wt) - imagE(N+1:2*N)*sin(wt),'bo', 'filled');
    %scatter(u(1:N)'+ realEM(1:N)*cos(wt) - imagEM(1:N)*sin(wt),u(N+1:2*N)'+ realEM(N+1:2*N)*cos(wt) - imagEM(N+1:2*N)*sin(wt),'bo', 'filled');
    hold off
    
    %scatter(u(1:N)'+stretch*mg(1:N,mode1).*cos(ph(1:N,mode1)+wt),u(N+1:2*N)'+stretch*mg(N+1:2*N,mode1).*cos(ph(N+1:2*N,mode1)+wt),'go', 'filled');
    hold off
    axis([-1.2*max(u0),1.2*max(u0),-1.2*max(u0),1.2*max(u0)])
    %title(['Magnetron Mode: ' num2str(mode1) ', Frequency = ' num2str(Dp(mode1)*wz/2/pi/1000) ' kHz, Rotation Frequency = 64 kHz, \delta = 0.0036' ])
    %title(['Cyclotron Mode: ' num2str(mode1-N) ', Frequency = ' num2str(Dp(mode1)*wz/2/pi/1000) ' kHz, Rotation Frequency = 64 kHz, \delta = 0.0036' ])
    %quiver(u(1:N)',u(N+1:2*N)',stretch*mg(1:N,mode).*cos(ph(1:N,mode)+wt),...
    %    stretch*mg(N+1:2*N,mode).*cos(ph(N+1:2*N,mode)+wt),0,'k','LineWidth',1)
        
      %  title(['Cyclotron Mode: ' num2str(mode1)])
    %writeVideo(vwr,getframe(gcf));
    %pause(.000001)
end
end
%close(vwr)
%pause(3)
%clf


%% Excite Planar Modes 2014_6_13

%date = '2014_6_23_PlanarExcitation';
%date = '2014_6_30_PlanarExcitation';
%date = '2014_7_30_PlanarExcitation_2';
%date = '2014_8_8_PlanarExcitation_2';
%date = '2014_8_13_PlanarExcitation';
%date = '2014_9_11_PlanarExcitation';
date = '2014_9_25_PlanarExcitation';
mkdir(['D:\PenningSimulationData\' date])
N = 19;  % Number of ions
f = 54.5;
Vwall = -80;
setTrapParameters(f,Vwall,N);% Set global variables
%N = 127;  % Number of ions
%setTrapParameters(44,-80,N);% Set global variables
global M l0 q m t0
Mass = diag([M M]);

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium

%Ep = Em;
%Dp = Dm;
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction
%Ep(2*N+1:end,:) = 1i*Ep(1:2*N,:).*repmat(Dp,1,2*N)'; % enforce eigenvectors to obey time derivative
%[Ep,Dp,st,V,A] = normalModes_OLD(u0,0); % Planar Modes with correction

for pMode = 1:2*N;
%for pMode = [1];
realE = real(Ep(1:2*N,pMode));
imagE = imag(Ep(1:2*N,pMode));


realEdot = real(Ep(2*N+1:end,pMode));
imagEdot = imag(Ep(2*N+1:end,pMode));
%unew = [u(1:N)+ realE(1:N)'*cos(wt) - imagE(1:N)'*sin(wt), u(N+1:2*N)+ realE(N+1:2*N)'*cos(wt) - imagE(N+1:2*N)'*sin(wt)];

scaling = 1e-0;
unew = u0 + scaling*realE'; 

%unew = [u(1:N)+ .1*realE(1:N)', u(N+1:2*N)+ 0.1*realE(N+1:2*N)']; % just excite positions?
%unew = u;

unew = unew*l0;
M = zeros(8,N);
M(1,:) = unew(1:N);
M(2,:) = unew(N+1:end);
M(3,:) = zeros(1,N);
%M(4,:) = -scaling*Dp(pMode)*imagE(1:N)';
%M(5,:) = -scaling*Dp(pMode)*imagE(N+1:end)';
M(4,:) = (l0/t0)*scaling*realEdot(1:N)';
M(5,:) = (l0/t0)*scaling*realEdot(N+1:end)';
M(7,:) = q*ones(1,N);
M(8,:) = m*ones(1,N);

filename = ['D:\PenningSimulationData\' date '\' num2str(N) '_' num2str(f) '_' num2str(Vwall) '_Pmode' num2str(pMode) '.dat'];
dlmwrite(filename,M,' ');
end

%% Load Planar Excitations and plot PSD (2014-6-17)

%FileLocation = 'D:\PenningSimulationData\2014_4_15_SmallCrystalModes\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_19_PlanarExcitation_3\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_20_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_23_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_31_PlanarExcitation\';  

%FileLocation = 'D:\PenningSimulationData\2014_8_6_PlanarExcitation\';  

%FileLocation = 'D:\PenningSimulationData\2014_6_30_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_9_11_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_9_19_PlanarExcitation\';  
FileLocation = 'D:\PenningSimulationData\2014_10_06_PlanarTemp\';  

setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);

freq = (0:0.5/(params(6)*params(7))/(params(5)/2):0.5/(params(6)*params(7)));
freq = 1.0e-6*freq(1:end-1); % chop of last one, because of Matlab ranges...

v = zeros(params(5),2*params(1));

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction
%Ep = Em;
%Dp = imag(Dm);
cmp = colormap(hsv(2*params(1)));
%cmp = lines;
figure
%for mode1 = 1:N
%for mode1 = N+1:2*N
%for mode1 = 1:2*N
for mode1 = [-1]    
[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode1);
%[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);
    
    
    for i = 0:params(5)-1
        %v(i+1,:) = [vxs(i+1,:) vys(i+1,:)]; % spin down velocities
        v(i+1,:) = spin_down(us(i+1,:),[vxs(i+1,:) vys(i+1,:)],params(2)); % spin down velocities
        us(i+1,:) = rotate(us(i+1,:),-thetas(i+1));           
    end
    % Calculate PSD for Planar Motion
    motion = us - repmat(u0,params(5),1); % subtract off equilibrium positions
    %motion = v;
    spectra = abs(fft(motion)).^2;
    Ppsd = sum(spectra, 2);
    Ppsd = Ppsd(1:(length(Ppsd)/2))+Ppsd(end:-1:length(Ppsd)/2+1);

    semilogy(freq,Ppsd,'Color',cmp(mode1,:))
    %semilogy(freq,Ppsd)
    xlabel('Frequency MHz','FontSize',24)
    ylabel('Planar PSD','FontSize',24)
    set(gca,'FontSize',24)
    
    disp('Done Calculating PSDs')
    


    
    hold on
    pause(.01)
end
   
%     
for i=1:2*params(1)
    plot([1e-6*wz/(2*pi)*Dp(i) 1e-6*wz/(2*pi)*Dp(i)],[1e-5,1e15],'k')
end

%%

mode = 2;
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction

        
InitialEnergy = real((Ep(1:2*N,mode))'*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*Ep(1:2*N,mode)'*V*Ep(1:2*N,mode));
EnergyMode=[];
for k = 1%params(5)   % I was wrong, this other linearization is not conventionally orthonormal
         
   %EnergyMode2(k) = 2*(us(k,:)*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*us(k,:)*V*Ep(1:2*N,mode))/Energy2;
   EnergyMode(k) = 2*(us(k,:)*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*us(k,:)*V*Ep(1:2*N,mode))/InitialEnergy;
   
end
%%
kB = 1.381e-23;
%mode = 8;
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction
%InitialEnergy = real((Ep(1:2*N,mode))'*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*Ep(1:2*N,mode)'*V*Ep(1:2*N,mode));
%EnergyMode = [];

%norm_coords_planar_test = 0.85*norm_coords_planar;
EnergyMode = zeros(params(5),2*params(1));
%EnergyMode2 = zeros(params(5),2*params(1));
for j=1:2*N
%for k = 1:params(5)   % I was wrong, this other linearization is not conventionally orthonormal
         
  
   %EnergyMode(k) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(k,mode))^2*(real((Ep(1:2*N,mode))'*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*(Ep(1:2*N,mode))'*V*(Ep(1:2*N,mode))))/InitialEnergy;
   %EnergyMode(k) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(k,mode))^2*(real((Ep(1:2*N,mode))'*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*(Ep(1:2*N,mode))'*V*(Ep(1:2*N,mode))));
   
  %!!! EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));%/real(EnergyP(j,j));
   EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,2*j-1)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));%/real(EnergyP(j,j));
   
   %EnergyMode2(:,j) = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)).^2*(real(( real(Ep(2*N+1:end,j)))'*real(Ep(2*N+1:end,j))  -0.5*real(Ep(1:2*N,j))'*V*real(Ep(1:2*N,j))));
   %EnergyMode2(:,j) = (1/4)*0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)).^2*(real(( 2*(Ep(2*N+1:end,j)))'*(Ep(2*N+1:end,j))  - (Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
   
   %EnergyMode2(:,j) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
   %EnergyMode2(:,j) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar_test(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
   %EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(:,j)).^2*((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j)));
   %EnergyMode2(:,j) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(:,j)).^2*(real((Ep(2*N+1:end,j))'*(Ep(2*N+1:end,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
    
end

%% Load Planar Excitations project into planar basis (2014-6-20)

%FileLocation = 'D:\PenningSimulationData\2014_6_19_PlanarExcitation_3\';  
FileLocation = 'D:\PenningSimulationData\2014_6_23_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_6_30_PlanarExcitation\';  

setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium
[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes with correction

kB = 1.381e-23;
v = zeros(params(5),2*params(1));
V_0 = PotentialPenning( u0);
%for mode1 = 1:N
%for mode1 = N+1:2*N
%for mode = 1:2*N
for mode = 8
   
    [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode);
    
    for i = 0:params(5)-1
    %for i = 0:10000-1
        %v(i+1,:) = [vxs(i+1,:) vys(i+1,:)]; % spin down velocities
        v(i+1,:) = spin_down(us(i+1,:),[vxs(i+1,:) vys(i+1,:)],params(2)); % spin down velocities
        us(i+1,:) = rotate(us(i+1,:),-thetas(i+1));      
        PPE(i+1) = 0.5*m*wz^2*l0^2*(PotentialPenning( us(i+1,:)) - V_0)/kB;
    end
    
    % Calculate PSD for Planar Motion
    motion = us - repmat(u0,params(5),1); % subtract off equilibrium positions
        
    %Project Planar Motion into Normal Modes
    norm_coords_planar = zeros(params(5),2*params(1)); 
    %norm_coords_planar = zeros(params(5),4*params(1)); 
   
    %InitialEnergy = real((Ep(1:2*N,mode))'*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*Ep(1:2*N,mode)'*V*Ep(1:2*N,mode));
    
    for k = 1:params(5)
    %for k = 1:10000
        %norm_coords_planar(k,:) = 2*(us(k,:)*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*us(k,:)*V*Ep(1:2*N,mode))/InitialEnergy;
        %norm_coords_planar(k,:) = 2*(motion(k,:)*(Dp(mode)*Dp(mode))*(Ep(1:2*N,mode))  -0.5*us(k,:)*V*Ep(1:2*N,mode))/InitialEnergy;
        temp = (inv(Ep)* [motion(k,:) t0*v(k,:)/l0]')';
        %norm_coords_planar(k,:) =  (inv(Ep)* [motion(k,:) t0*v(k,:)/l0]')';
        norm_coords_planar(k,:) =  temp(1:2:end) + temp(2:2:end) + (temp(1:2:end) - temp(2:2:end)) ;
        PPE2(k) = -0.5*m*wz^2*l0^2*motion(k,:)*V*motion(k,:)'/kB;

    end

   
    %save([FileLocation 'planarMode_' num2str(mode) '.mat'],'norm_coords_planar')
    
end
   %% Load Planar Laser Cooling (2014-10-16) simple

%FileLocation = 'D:\PenningSimulationData\2014_5_11_SmallCrystalModes\'; 
%FileLocation = 'D:\PenningSimulationData\2014_6_22_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_7_15_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_6_23_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_29_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_30_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_13_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_6_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_10_11_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_10_28_PlanarTemp\'; 
FileLocation = 'D:\PenningSimulationData\2014_11_05_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_25_PlanarExcitation\';  

setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
kB = 1.381e-23;

[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
[u0 V_0] = findEquilibrium(u0);    % Solve for equilibrium
%[u0 sigma] = normal_traj(us,u0);   %shifted u0
[Ea,Da] = normalModes(u0,1); % Axial Modes
[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes

norm_coords = zeros(params(5),params(1)); 
norm_coords_planar = zeros(params(5),4*params(1)); 
v = zeros(params(5),2*params(1));
EnergyMode = zeros(params(5),2*params(1));

cmp = colormap(hsv(N));
%V_0 = PotentialPenning( rotate(us(1,:),-thetas(1)));
%V_0 = 0;
%Project Axial and Planar Motion in to Normal Modes
for i = 1:params(5)
    
    vlab(i,:) = [vxs(i,:) vys(i,:)]; % lab velocities
    v(i,:) = spin_down(us(i,:),[vxs(i,:) vys(i,:)],params(2)); % spin down velocities
    us(i,:) = rotate(us(i,:),-thetas(i));  
    PPE(i) = 0.5*m*wz^2*l0^2*( PotentialPenning(us(i,:)) - V_0 )/kB;
    norm_coords(i,:) = modeBasis(zs(i,:),Ea); 
    motion = us(i,:) - u0; % subtract off equilibrium positions  
    norm_coords_planar(i,:) = (inv(Ep)* [motion t0*v(i,:)/l0]')';
    
end
PKE = 0.5*m*sum(v.^2,2)/kB;
for k = 1:2*N
    EnergyMode(:,k) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2);  
   % EnergyMode(:,k) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2)*( (real((Ep(1:2*N,2*k-1))'*(Dp(2*k)*Dp(2*k))*(Ep(1:2*N,2*k-1))  -0.5*Ep(1:2*N,2*k-1)'*V*Ep(1:2*N,2*k-1))) + ...
    %                    (real((Ep(1:2*N,2*k))'*(Dp(2*k)*Dp(2*k))*(Ep(1:2*N,2*k))  -0.5*Ep(1:2*N,2*k)'*V*Ep(1:2*N,2*k))));
end
save([FileLocation 'axialModeDecomposition.mat'],'norm_coords')
save([FileLocation 'planarModeDecomposition.mat'],'norm_coords_planar','EnergyMode')
save([FileLocation 'DataVars.mat'],'PKE','PPE')

%% Load Planar Laser Cooling (2014-6-24)

%FileLocation = 'D:\PenningSimulationData\2014_5_11_SmallCrystalModes\'; 
%FileLocation = 'D:\PenningSimulationData\2014_6_22_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_7_15_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_6_23_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_29_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_7_30_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_8_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_13_PlanarExcitation\';  
%FileLocation = 'D:\PenningSimulationData\2014_8_6_PlanarExcitation_2\';  
%FileLocation = 'D:\PenningSimulationData\2014_10_11_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarTemp\'; 
FileLocation = 'D:\PenningSimulationData\2014_10_28_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_25_PlanarExcitation\';  

setTrapParameters(0,0,0);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
kB = 1.381e-23;

[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);
%[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,2);
%for mode=1:2*N
%for mode=2
    %mode
%[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode);


u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
[u0 V_0] = findEquilibrium(u0);    % Solve for equilibrium

PPE1 = zeros(1,params(5));
for i = 1:params(5)
   
    %vlab(i,:) = [vxs(i,:) vys(i,:)]; % lab velocities
    v(i,:) = spin_down(us(i,:),[vxs(i,:) vys(i,:)],params(2)); % spin down velocities
    us(i,:) = rotate(us(i,:),-thetas(i));  
    PPE1(i) = 0.5*m*wz^2*l0^2*( PotentialPenning(us(i,:)) - V_0 )/kB;  
    
end
[u0 sigma] = normal_traj(us,u0);   %shifted u0
[Ea,Da] = normalModes(u0,1); % Axial Modes

[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes

norm_coords = zeros(params(5),params(1)); 
norm_coords_planar = zeros(params(5),4*params(1)); 
v = zeros(params(5),2*params(1));
EnergyMode = zeros(params(5),2*params(1));

cmp = colormap(hsv(N));
%V_0 = PotentialPenning( rotate(us(1,:),-thetas(1)));
%V_0 = 0;
%Project Axial and Planar Motion in to Normal Modes
PPE2 = zeros(1,params(5));
V_2 = PotentialPenning(u0);
for i = 1:params(5)
    
    %vlab(i,:) = [vxs(i,:) vys(i,:)]; % lab velocities
    %v(i,:) = spin_down(us(i,:),[vxs(i,:) vys(i,:)],params(2)); % spin down velocities
    %us(i,:) = rotate(us(i,:),-thetas(i));  
    PPE2(i) = 0.5*m*wz^2*l0^2*( PotentialPenning(us(i,:)) - V_2 )/kB;
    %norm_coords(i,:) = modeBasis(zs(i,:),Ea); 
    %motion = us(i,:) - u0; % subtract off equilibrium positions
    
    %Ep(:,1:2) = -Ep(:,1:2);
    %Ep(:,1:2) = -real(Ep(:,1:2));
    %temp = (inv(Ep)* [motion t0*v(i,:)/l0]')';
    %norm_coords_planar(i,1) = dot(motion,Ep(1:2*N,1));
    %test(i) = dot(motion,-Ep(1:2*N,1));
    
    %norm_coords_planar(i,:) =  temp(1:2:end) + temp(2:2:end) + (temp(1:2:end) - temp(2:2:end)); %store real and imaginary components
   
    %norm_coords_planar(i,:) = (inv(Ep)* [motion t0*v(i,:)/l0]')';
    
end
% PKE = 0.5*m*sum(v.^2,2)/kB;
% [Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes
% for k = 1:2*N
%     j = k;
%  %   EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,2*j-1)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));%/real(EnergyP(j,j));
%    EnergyMode(:,k) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2);
%    % EnergyMode(:,k) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*k-1)).^2+abs(norm_coords_planar(:,2*k)).^2)*( (real((Ep(1:2*N,2*k-1))'*(Dp(2*k)*Dp(2*k))*(Ep(1:2*N,2*k-1))  -0.5*Ep(1:2*N,2*k-1)'*V*Ep(1:2*N,2*k-1))) + ...
%    %         (real((Ep(1:2*N,2*k))'*(Dp(2*k)*Dp(2*k))*(Ep(1:2*N,2*k))  -0.5*Ep(1:2*N,2*k)'*V*Ep(1:2*N,2*k))));
% end

% load([FileLocation 'planarMode_' num2str(mode) '.mat'],'norm_coords_planar')
% load([FileLocation 'DataVars' num2str(mode) '.mat'],'PKE','PPE')
% EnergyMode = zeros(params(5),4*params(1));
% for j=1:4*N
%     %EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
%     EnergyMode(:,j) = (real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
% end
% time = (0:params(5)-1)*(params(7)*params(6));

% figure(10)
% %semilogx(time,real(norm_coords_planar(:,1)),'Color',cmp(mode,:))
% semilogx(time,real(EnergyMode(:,1)),'Color',cmp(mode,:))
% hold on
% figure(11)
% loglog(time,real(EnergyMode(:,mode)),'Color',cmp(mode,:))
% hold on
% figure(12)
% loglog(time,PKE,'Color',cmp(mode,:))
% hold on
%save([FileLocation 'axialMode_' num2str(mode) '.mat'],'norm_coords')
%save([FileLocation 'planarMode_' num2str(mode) '.mat'],'norm_coords_planar','EnergyMode')
%save([FileLocation 'DataVars' num2str(mode) '.mat'],'PKE','PPE')
%end
%save([FileLocation 'axialModeDecomposition.mat'],'norm_coords')
%save([FileLocation 'planarModeDecomposition.mat'],'norm_coords_planar','EnergyMode')
%save([FileLocation 'DataVars.mat'],'PKE','PPE')

%% this doesn't work right
[Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes
ModeEnergy = zeros(params(5),2*params(1));
for j = 1:2*params(1)
    %ModeEnergy(:,j) = 0.5*m*(wz*Dp(j))^2*l0^2*(real(norm_coords_planar(:,j)))'.^2/(kB)/real(EnergyP(j,j));
    ModeEnergy(:,j) = 0.5*m*(wz*Dp(j))^2*l0^2*(abs(norm_coords_planar(:,j)))'.^2/(kB);
    %ModeEnergy(:,j) = 0.5*m*(wz)^2*l0^2*(real(norm_coords_planar(:,j)))'.^2/(kB)/real(EnergyP(j,j));
    %ModeEnergy(:,j) = 0.5*m*(wz)^2*l0^2*(real(norm_coords_planar(:,j)))'.^2/kB;
    %ModeEnergy(:,j) = 0.5*m*(wz*Dp(j))^2*(real(norm_coords_planar(:,j))'.^2)/(kB);
end
%%

AxialModeEnergy = zeros(params(5),params(1));
for j = 1:params(1)
    %AxialModeEnergy(:,j) = 0.5*m*(wz*D(j))^2.*norm_coords(:,j)'.^2/(kB/2);
    AxialModeEnergy(:,j) = 0.5*m*(wz*Da(j))^2.*norm_coords(:,j)'.^2/kB;
end

%% 

total = sum(EnergyMode1,2);
endiff = total - (PPE' + PKE);


%% Non linear effects planar modes 8/13/14

setTrapParameters(64,-80,7);
%setTrapParameters(54.5,400,19);
%setTrapParameters(64,-1600,7);
%setTrapParameters(44,-80,127);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
[u0, V_0] = findEquilibrium(u0);    % Solve for equilibrium
%u0 = rotate(u0,pi/2);       % get same major axis as dominics code

%setTrapParameters(54.5,-400,19);
%setTrapParameters(44,-80,127);
%u0 = rotate(u0,pi/2);       % get same major axis as dominics code

[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes
hbar = 1.05457173e-34;
kB = 1.381e-23;
Tdoppler = hbar*2*pi*19e6/2/kB;
gz = (Dp>0);   % get positive values
Dp = Dp(gz);

phi = linspace(0,pi/2,10);

fig1 = figure;
cmp = colormap(hsv(N));
fig2 = figure;
fig3 = figure;
cmp = colormap(hsv(N));
xlog = logspace(-10,10,1000);

%for mode=13:14
%for mode=N+1:2*N
    LinearEnergy = [];
    IonEnergy = [];
for mode=1:N
%for mode=N
    index = 1;

    for scale = xlog
        
        norm_coords_planar = zeros(1,4*N);
        norm_coords_planar(2*mode-1:2*mode) = -scale*[1 1];
        %norm_coords_planar(1:2) = (1/82.648774640897884)*[1 1];
        
        %norm_coords_planar(2*mode-1:2*mode) = scale*[(1+1i)/sqrt(2) (1+1i)/sqrt(2)];
        %norm_coords_planar(2*mode-1:2*mode) = scale*[exp(1i*phi) exp(1i*phi)];
        %EnergyMode = zeros(1,2*params(1));
        %for j=1:2*N
        %    EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(2*j-1)).^2+abs(norm_coords_planar(2*j)).^2);
        %end
        
        EnergyMode = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(2*mode-1)).^2+abs(norm_coords_planar(2*mode)).^2);
        %total = sum(EnergyMode(:,1:end),2);
        
        ions = (Ep*norm_coords_planar.');
        positions = ions(1:2*N)' + u0;
        velocities = l0*ions(2*N+1:end)/t0;
        %PPE = 0.5*m*wz^2*l0^2*( PotentialPenning( rotate(positions,-pi/2)) - V_0 )/kB;
        PPE = 0.5*m*wz^2*l0^2*( PotentialPenning( positions) - V_0 )/kB;
        PKE = 0.5*m*sum(velocities.^2)/kB;
        IonEnergy(index) = PPE+PKE;
        LinearEnergy(index) = EnergyMode;
        index = index + 1;
    end
    %index
    %figure
    %loglog(LinearEnergy,LinearEnergy,'k')
    %hold on
    %loglog(LinearEnergy,IonEnergy,'m')
    
    EnergyGap = abs(LinearEnergy - IonEnergy)/Tdoppler;
    %EnergyGap = abs(LinearEnergy./IonEnergy);
    %semilogy(EnergyGap,'Color',cmp(mode,:))
    %loglog(xlog,EnergyGap,'Color',cmp(mod(mode,N+1)+1,:))
    %loglog(fig1,LinearEnergy/Tdoppler,EnergyGap/Tdoppler,'Color',cmp(mode-N,:))
    figure(fig1)
    loglog(LinearEnergy/Tdoppler,EnergyGap,'Color',cmp(mode,:))
    %loglog(IonEnergy/Tdoppler,EnergyGap/Tdoppler,'Color',cmp(mode,:))
    hold on
    
    
    %loglog(fig1,LinearEnergy/Tdoppler,(hbar*Dp(mode)*wz/2/pi/kB)*ones(size(EnergyGap))/Tdoppler,'--','Color',cmp(mode-N,:))
   
    loglog(LinearEnergy/Tdoppler,(hbar*Dp(mode)*wz/2/pi/kB)*ones(size(EnergyGap))/Tdoppler,'--','Color',cmp(mode,:))
    %loglog([modeTemp(mode) modeTemp(mode)]/Tdoppler,[10^-20 10^20],'-.','Color',cmp(mode,:))
    %loglog(IonEnergy/Tdoppler,(hbar*Dp(mode)*wz/2/pi/kB)*ones(size(EnergyGap))/Tdoppler,'--','Color',cmp(mode,:))
    loglog([Tdoppler Tdoppler]/Tdoppler,[10^-20 10^20],'k:')
    loglog([10^-20 10^20],[Tdoppler Tdoppler]/Tdoppler,'k:')

    %axis([10^-20,10^20,10^-20,10^20])
    axis([10^-15,10^15,10^-15,10^15])
    
    xlabel('Linear Theory Energy/T_{Doppler}','FontSize',24)
    ylabel('Absolute Energy Difference between Linear and Exact/T_{Doppler}','FontSize',24)
    
    %xlabel('Linear Theory Energy (Kelvin)','FontSize',24)
    %ylabel('Absolute Energy Difference between Linear and Exact (Kelvin)','FontSize',24)
    %title(['Magnetron Branch Nonlinearities ' num2str(Vw)] , 'FontSize',24)
    title(['Magnetron Branch Nonlinearities '] , 'FontSize',24)
    %title('Cyclotron Branch Nonlinearities', 'FontSize',24)
    clb = colorbar;
    set(clb,'YTick',1:N)
    set(gca,'FontSize',20)
    
    figure(fig3)
    EnergyGap = abs(LinearEnergy./IonEnergy);
    loglog(LinearEnergy/Tdoppler,EnergyGap,'Color',cmp(mode,:))
    hold on
    loglog([Tdoppler Tdoppler]/Tdoppler,[10^-20 10^20],'k:')
    axis([10^-15,10^15,10^-15,10^15])
    xlabel('Linear Theory Energy/T_{Doppler}','FontSize',24)
    ylabel('Ratio of Linear Theory to Exact','FontSize',24)
    title(['Magnetron Branch Nonlinearities '] , 'FontSize',24)
    clb = colorbar;
    set(clb,'YTick',1:N)
    set(gca,'FontSize',20)
    
    figure(fig2)
    EnergyGap = abs(LinearEnergy - IonEnergy)/Tdoppler;
    %if mode == 2
    %    continue
    %end
    plot(diff(log(EnergyGap/Tdoppler))./diff(log(LinearEnergy/Tdoppler)),'Color',cmp(mode,:));
    hold on
    ylabel('Slope of LogLog plot','FontSize',24)
    
    title(['Magnetron Branch Nonlinearities Slopes'] , 'FontSize',24)
        
end
%% Calculate Planar Temperatures
%FileLocation = 'D:\PenningSimulationData\2014_7_15_PlanarTemp\';
%FileLocation = 'D:\PenningSimulationData\2014_7_17_PlanarTemp\'; 
%FileLocation = 'D:\PenningSimulationData\2014_9_2_PlanarTemp\';
FileLocation = 'D:\PenningSimulationData\2014_10_08_PlanarTemp\';
%FileLocation = 'D:\PenningSimulationData\2014_9_30_PlanarTemp\';


load([FileLocation 'DataVars.mat'],'PKE','PPE')
load([FileLocation 'planarModeDecomposition.mat'],'norm_coords_planar')
setTrapParameters(0,0,0)
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(params(2),-params(3)/G,params(1));
thetas = dlmread([FileLocation 'thetas.dat']);
times = linspace(0,params(5)*params(7)*params(6),params(5));

kB = 1.381e-23;
hbar = 1.05457173e-34;
Tdoppler = hbar*2*pi*19e6/2/kB;

%[us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);

%u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
%u0 = rotate(u0,pi/2);       % get same major axis as dominics code
%u0 = findEquilibrium(u0);    % Solve for equilibrium
%[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes

EnergyMode = zeros(params(5),2*params(1));
for j=1:2*N
    EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*(abs(norm_coords_planar(:,2*j-1)).^2+abs(norm_coords_planar(:,2*j)).^2);
end

totalEmode = sum(EnergyMode(:,2:end),2); % don't include 1st mode because it's wrong
modeTemp = [];
cutoff = 10; %remove first tenth to calculate temperature
for j = 1:2*N
    modeTemp(j) = mean(EnergyMode(length(EnergyMode)/cutoff:end,j));
    %modeTemp(j) = 0.5*m*(wz*D(j))^2*mean(norm_coords(1000:end,j)'.^2)/(kB/2);
end
totalE = PPE+PKE';
figure



%scatter(1:length(modeTemp),modeTemp/Tdoppler)
scatter(2:length(modeTemp),modeTemp(2:end)/Tdoppler,'filled')
hold on
plot([0 2*N], mean(totalE(length(EnergyMode)/cutoff:end)/(2*N))*[1 1]/Tdoppler,'m')
plot([0 2*N], mean(totalEmode(length(EnergyMode)/cutoff:end)/(2*N))*[1 1]/Tdoppler,'g')
title(['Planar Temperature, 19 ions, dt = 5e-10, Vwall = -400, \omega = 54.5 '] , 'FontSize',24)
% scatter(2:length(modeTemp),modeTemp(2:end)/1e-3,'filled')
% hold on
% plot([0 2*N], [Tdoppler Tdoppler]/1e-3,'k:')
% plot([0 2*N], mean(totalE(length(EnergyMode)/10:end)/(2*N))*[1 1]/1e-3,'m')
% plot([0 2*N], mean(totalEmode(length(EnergyMode)/10:end)/(2*N))*[1 1]/1e-3,'g')

set(gca,'FontSize',24)
xlabel('Mode Number','FontSize',24)
%ylabel('Temperature (mK)','FontSize',24)
ylabel('Temperature/T_{Doppler}','FontSize',24)
%%
index = 1;
for Vwall = linspace(-1000,-2000,100)
    setTrapParameters(64,Vwall,7);
    global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
    u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
    u0 = rotate(u0,pi/2);       % get same major axis as dominics code
    [u0 V_0] = findEquilibrium(u0);    % Solve for equilibrium
    [Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes
    test(index) = Dp(1);
    index = index + 1;
end

%% Non linear effects axial modes

%setTrapParameters(64,-80,7);
setTrapParameters(54.5,-80,19);
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params t0
u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
[u0 V_0] = findEquilibrium(u0);    % Solve for equilibrium
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
[Ea,Da] = normalModes(u0,1); 

fig1 = figure;
cmp = colormap(hsv(N));
fig2 = figure;
    LinearEnergy = [];
    IonEnergy = [];
for mode=1:N
    index = 1;
    xlog = logspace(-30,20,10000);
    
    for scale = xlog
        norm_coords = zeros(1,N);
        norm_coords(mode) = scale;
        z = Ea*norm_coords'; % scale 
        
        EnergyMode = 0.5*m*(wz*Da(mode))^2*l0^2*(norm_coords(mode)'.^2)/kB;
        r0 = pairFullDistance(u0,zeros(1,N)); % Calculate pairwise distances
        ZPE0 = 1./r0;                                  % Find Coulomb potential between ions
        ZPE0(isnan(ZPE0) | isinf(ZPE0)) = 0;            % Self potential = 0
        ZPE0 =  0.5*m*wz^2*l0^2*(0.5*sum(sum(ZPE0)))/kB;        % Total potential energy
        
        r = pairFullDistance(u0,z); % Calculate pairwise distances
	    ZPE = 1./r;                                  % Find Coulomb potential between ions
        ZPE(isnan(ZPE) | isinf(ZPE)) = 0;            % Self potential = 0
        ZPE =  0.5*m*wz^2*l0^2*(sum(z.^2) + 0.5*sum(sum(ZPE)))/kB - ZPE0;        % Total potential energy
        %ZPE =  0.5*m*wz^2*l0^2/kB*(sum(z.^2));        % Total potential energy
       % ZPE =  0.5*m*wz^2*l0^2*(0.5*sum(sum(ZPE))) %-ZPE0;        % Total potential energy
        
        IonEnergy(index) = ZPE;
        LinearEnergy(index) = EnergyMode;
        index = index + 1;
    end
%     loglog(xlog,LinearEnergy,'b')
%     hold on
%     loglog(xlog,IonEnergy,'g')

    EnergyGap = abs(LinearEnergy - IonEnergy);
    figure(fig1)
    loglog(LinearEnergy/Tdoppler,EnergyGap/Tdoppler,'Color',cmp(mode,:))
    hold on
    loglog(LinearEnergy/Tdoppler,(hbar*Da(mode)*wz/2/pi/kB)*ones(size(EnergyGap))/Tdoppler,'--','Color',cmp(mode,:))
    loglog([Tdoppler Tdoppler]/Tdoppler,[10^-30 10^30],'k:')
    loglog([10^-30 10^30],[Tdoppler Tdoppler]/Tdoppler,'k:')
    axis([10^-30,10^30,10^-30,10^30])
    clb = colorbar;
    set(clb,'YTick',1:N)
    set(gca,'FontSize',20)
    xlabel('Linear Theory Energy/T_{Doppler}','FontSize',24)
    ylabel('Absolute Energy Difference between Linear and Exact/T_{Doppler}','FontSize',24)
    title(['Axial Branch Nonlinearities '] , 'FontSize',24)
    
    figure(fig2)
    if mode == 7
        
        continue
    end
    plot(diff(log(EnergyGap/Tdoppler))./diff(log(LinearEnergy/Tdoppler)),'Color',cmp(mode,:));
    hold on
    ylabel('Slope of LogLog plot','FontSize',24)
    
    title(['Axial Branch Nonlinearities Slopes'] , 'FontSize',24)
end
%%
setTrapParameters(54.5,-400,19);
for phi= [0 pi/2 pi 3*pi/2]
    u0 = generateLattice(N,1);
    u0 = rotate(u0,phi);   
    [u0 V_0] = findEquilibrium(u0);
    [phi/pi V_0]
end
%% 

stepssizes = [5e-8, 5e-9, 5e-10, 5e-11];
colors = ['m'; 'r'; 'b'; 'g'];
files ={'2014_8_5_PlanarExcitation'; '2014_7_31_PlanarExcitation'; '2014_8_1_PlanarExcitation'; '2014_8_4_PlanarExcitation_2'};

for i = 3%1:4
    FileLocation = ['D:\PenningSimulationData\' files{i} '\'];
    load([FileLocation 'planarModeDecomposition.mat'],'norm_coords_planar')
    load([FileLocation 'DataVars.mat'],'PKE','PPE')
    time = (0:1e6-1)*(stepssizes(i)*100);
    [Ep,Dp,st,V,A] = normalModes(u0,0); % Planar Modes with correction
    EnergyMode = zeros(1e6,2*params(1));
    for j=1:2*N
        EnergyMode(:,j) = 2*0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)/2).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
    end
    j=1;
    EnergyMode(:,j) = 0.5*m*wz^2/kB*l0^2*real(norm_coords_planar(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));
    total = sum(EnergyMode(:,2:end),2);
    %EnergyMode1 = 0.5*m*(wz*Dp(1))^2*l0^2*(abs(norm_coords_planar(:,1)))'.^2/(kB);
    j=2;
%    EnergyMode1 = 0.5*m*wz^2/kB*l0^2*abs(norm_coords_planar(:,j)).^2*(real((Ep(1:2*N,j))'*(Dp(j)*Dp(j))*(Ep(1:2*N,j))  -0.5*(Ep(1:2*N,j))'*V*(Ep(1:2*N,j))));

    %figure(1)
    %plot(time,real(norm_coords_planar(:,1)),colors(i))
    %hold on
    
    figure(7)
    %semilogx(time,real(EnergyMode1),colors(i))
    %loglog(time,real(EnergyMode1),colors(i))
    loglog(time,real(total),colors(i))
    %semilogx(time,real(norm_coords_planar(:,1)),colors(i))
    hold on
    
    figure(8)
    %loglog(time,PKE,colors(i))
    %hold on
    loglog(time,PKE'+PPE-0.5*m*wz^2*l0^2*V_0/kB,'b')
    hold on
    loglog(time,real(total),'k')
    loglog(time,EnergyMode(:,1),'y')
     
    figure(9)
    %plot(time,PPE-0.5*m*wz^2*l0^2*V_0/kB,colors(i))
    hold on
    
    %figure(4)
    %semilogy(time,PKE,colors(i))   
    %hold on
    
end


%% Find frequency of mode from data
Fs = 1/(100*5e-9);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(real(norm_coords_planar(:,1)),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(Y(1:NFFT/2+1))) 