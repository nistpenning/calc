% Quickly display snapshots of simulation to make a "movie"

cleanupObj = onCleanup(@CLEAN);

%FileLocation = 'C:\Users\ACKWinDesk\Documents\GitHub\ultracold-ions-notes\axialPhononSpectra\data\run4extra\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_19_EquilibriumTest\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_20_NormalModeExpansionTest_BigData\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_26_PlaneTransition\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_26_HighFrequencyTest\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_23_NormalModeExpansionTest_BigData_smallZ-8\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_5_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_9_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_27_NormalModeExpansionTest_inPlane\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_7_AxialLaserCooling_COM1e-5\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_25_AxialTemp_LaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_8_EnergyTest\';
%FileLocation = 'D:\PenningSimulationData\2014_3_28_SmallCrystalModes\';    
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_5_NonEquilibriumTest_Mode45_zvels\';
FileLocation = 'D:\PenningSimulationData\2014_3_27_AxialTemp_LaserCooling\';



%setTrapParameters(44,-80,127);
%global w % get rotation velocity
%dt = 5.0e-9;
%delta_theta = w*dt;
setTrapParameters(0,0,0);
thetas = dlmread([FileLocation 'thetas.dat']);
params = dlmread([FileLocation 'params.dat']);
global m wz l0 q V0 G 
setTrapParameters(params(2),-params(3)/G,params(1));
%filename =[FileLocation 'run4_0.dat'];
%filename =[FileLocation '0.dat'];
%M = dlmread(filename);
%u = convertPythonDataToMatlab(M);
%theta = findTheta(u);
%theta = 1.6228; %original run4
%theta = 0.9896;  % run4 extra
%theta = 0.5850 * pi;
%u = rotate(u,theta);            % rotate to that config
%u = rotate(u,thetas(1));            % rotate to that config

%figure(1)
%scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
%axis([-25,25,-25,25])
% figure(2)
% scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
% axis([-25,25,-25,25])
% pause(.01)

vidpath = 'OnlyAxialBeam';
%vidpath = '19ions';
vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
vwr.FrameRate = 30;
open(vwr);

%cmp = colormap(hsv(6));
% figure(2)
% scatter(u(1:end/2),u(end/2+1:end),20,cmp(1,:), 'filled');
% axis([-25,25,-25,25])

%for i = 0:4999
%for i = [0 999 1999 2999 3999 4999]
for i = 0:2000:params(5)-1
%for i = 0:1000:500000
%    theta = theta + delta_theta;
 %   filename =[FileLocation 'run4_' int2str(i) '.dat'];
    %filename =[FileLocation int2str(i) '_*.dat'];
    filename =[FileLocation int2str(i) '.dat'];
    M = dlmread(filename);
    [u z] = convertPythonDataToMatlab(M);
  %  figure(1)
   % scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
   % axis([-25,25,-25,25])
    
   % u = rotate(u,theta);            % rotate to that config
    u = rotate(u,-thetas(i+1));            % rotate to that config
    if i==0
        u0=u;
    end
   % figure(1)
     scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
     hold on
    % scatter(u0(1:end/2),u0(end/2+1:end),'go', 'filled');
     hold off
%     dists(i/1000+1) = ConfigurationDistance(u,u0);
%     hold on
%     scatter(u(end/2),u(end),'yo', 'filled');
%     
%     scatter(u(50),u(127+50),'co', 'filled');
%     scatter(u(1),u(end/2+1),'go', 'filled');
%     scatter(u(25),u(127+25),'mo', 'filled');
%     hold off

    %hold on
    %scatter(u(1:end/2),u(end/2+1:end),20,cmp(floor((i+1)/1000)+1,:), 'filled');
   % figure(1)
   % scatter(u(1:end/2),u(end/2+1:end),20,z/max(abs(z)), 'filled');
    axis([-25,25,-25,25])
    %axis([-9,9,-9,9])
    %title('Lab Frame, N = 127, delta = 0.0036, w = 44 kHz, ')
    title(['Rotating Frame, N = ' num2str(params(1)) ', \delta = ' num2str(params(3)) ', w = ' num2str(params(2)) ' kHz, i = ' num2str(i) ', t = ' num2str(i*params(6)*params(7)) ])
    xlabel('Scaled position')
   % figure(2)
    %scatter(i,z(100))
    %hold on
%     
%     if mod((i+1),1000)==0
%         figure(2)
%         hold on
%         scatter(u(1:end/2),u(end/2+1:end),1,cmp(floor((i+1)/1000),:), 'filled');
%         axis([-25,25,-25,25])
%     end
    pause(.00001)
   % pause(1)

    writeVideo(vwr,getframe(gcf));
   % clf

end
close(vwr);


% params(1) = 127;
% params(2) = 44;
% params(3) = 0.0036;
% params(4) = 0;
% params(5) = 5e5;
% params(6) = 100;
% params(7) = 5e-9;