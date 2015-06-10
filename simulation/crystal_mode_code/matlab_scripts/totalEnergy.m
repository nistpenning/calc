setTrapParameters(0,0,0);
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_9_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_16_AxialTemp_LaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_3_3_EnergyTest\';
%FileLocation = 'D:\PenningSimulationData\2014_3_25_AxialTemp_LaserCooling\';  
FileLocation = 'D:\PenningSimulationData\2014_3_27_AxialTemp_LaserCooling\';   

params = dlmread([FileLocation 'params.dat']);
thetas = dlmread([FileLocation 'thetas.dat']);

load([FileLocation 'axialModeDecomposition_eq.mat'],'norm_coords')
load([FileLocation 'axialVelModeDecomposition_eq.mat'],'norm_vels')

global m wz l0 q V0 G 
setTrapParameters(params(2),-params(3)/G,params(1));

filename = [FileLocation int2str(params(5)-1) '.dat']; 
M = dlmread(filename);
u = convertPythonDataToMatlab(M);
u = rotate(u,-thetas(params(5)-1)); 
u = findEquilibrium(u);                % original equilibrium lattice
[E,D,st,A] = normalModes(u,1);

%AxialPotentialIon = zeros(params(5),1);
%AxialKineticIon = zeros(params(5),1);
%AxialKineticModes = zeros(params(5),1);
%AxialPotentialModes = zeros(params(5),1);
%AxialTrapIon = zeros(params(5),1);
PlanarKineticIons = zeros(params(5),1);
PlanarKineticIons_rotframe = zeros(params(5),1);

for i = 1:params(5)

    filename =[FileLocation int2str(i-1) '.dat'];
    M = dlmread(filename);
    u = convertPythonDataToMatlab(M); % dimensionless positions
%     x = l0*u(1:end/2);
%     y = l0*u(end/2+1:end);
    z = M(3,:);
    vz = M(6,:);
    
     AxialTrapIon(i) = q*V0*sum(z.^2);
     AxialKineticIon(i) = 0.5*m*sum(vz.^2);
%     AxialPotentialIon(i) = 0.5*m*wz^2*(z*A*z');
%     AxialKineticModes(i) = 0.5*m*sum(norm_vels(i,:).^2);
%     AxialPotentialModes(i) = 0.5*m*sum( ((wz*D).^2).*(norm_coords(i,:)'.^2));
    
%      vx = M(4,:);
%      vy = M(5,:);
%      PlanarKineticIons(i) = 0.5*m*sum(vx.^2 + vy.^2);
%      v = spin_down(u,[vx vy],params(2));
%      PlanarKineticIons_rotframe(i) = 0.5*m*sum(v.^2);
     
end
%save([FileLocation 'Energies.mat'],'AxialPotentialIon','AxialKineticIon','AxialKineticModes','AxialPotentialModes')
%save([FileLocation 'PlanarEnergies.mat'],'PlanarKineticIons','PlanarKineticIons_rotframe')
% 
% params = dlmread([FileLocation 'params.dat']);
% thetas = dlmread([FileLocation 'thetas.dat']);
% 
% global m wz l0 q V0 ke B G Vw w
% binsize = 1000; 
% k = 1;
% setTrapParameters(params(2),-params(3)/G,params(1));
% filename =[FileLocation int2str(0) '.dat'];
% M = dlmread(filename);
% u0 = convertPythonDataToMatlab(M); % dimensionless positions
% minEnergy = 0.5*m*wz^2*l0^2*PotentialPenning(u0);
% E1Pot = zeros(binsize,127);
% E1KE = zeros(binsize,127);
% AxialPotEnergy = [];
% AxialKEnergy = [];
% 
% %for i = 1:binsize:params(5)
% for i = 1:1000
% 
%     ind = (k-1)*params(5)/binsize + i-1;
%     filename =[FileLocation int2str(ind) '.dat'];
%     %filename = [FileLocation int2str(i) '.dat']; 
%  %  i = (i-1)/binsize+1;
%     M = dlmread(filename);
%     u = convertPythonDataToMatlab(M); % dimensionless positions
%     x = l0*u(1:end/2);
%     y = l0*u(end/2+1:end);
%     
%     %v = [M(4,:) M(5,:)];              % get velocities  
%     vx = M(4,:);
%     vy = M(5,:);
%     PlanarKEnergy(i) = 0.5*m*sum(vx.^2 + vy.^2);
% %     theta = thetas(ind+1);
% %     V=1./pairDistance(l0*u);                 % Find Coulomb potential between ions
% %     V(isnan(V) | isinf(V))=0;                % Self potential = 0
% %     PlanarPotEnergy(i) =  0.5*ke*q^2*sum(sum(V)) ...
% %        -0.5*q*V0*sum(x.^2 + y.^2) + ...
% %         q*G*Vw*V0*sum((x.^2-y.^2)*cos(2*theta) - 2*x.*y*sin(-2*theta)) + ...
% %         0.5*q*B*sum(vx.*y - vy.*x);
% %     rotatingWall(i) = sum((x.^2-y.^2)*cos(2*theta) - 2*x.*y*sin(-2*theta));
% %     v = spin_down(u,v,params(2));
% %     u = rotate(u,-thetas(ind+1));     % rotate positions
% %     v = rotate(v,-thetas(ind+1));     % rotate velocities 
%      z = M(3,:);
%      vz = M(6,:);
%     % vz = M(4,:);
%       AxialPotEnergy(i) = q*V0*sum(z.^2);
%      AxialKEnergy(i) = 0.5*m*sum(vz.^2);
%      
%     % hist(vz)
%      %pause(.1)
% %     
% %     x = l0*u(1:end/2);
% %     y = l0*u(end/2+1:end);
% %     vx = v(1:end/2);
% %     vy = v(end/2+1:end);
% %     %AxialEnergy(i) = sum(ionAxialEnergy(z,vz));
% %    % E1Pot(i,:) = q*V0*z.^2;
% %     %E1KE(i,:) = 0.5*m*vz.^2;
% %     AxialPotEnergy(i) = q*V0*sum(z.^2);
% %     AxialKEnergy(i) = 0.5*m*sum(vz.^2);
% %      PlanarPotEnergy_rot(i) = 0.5*m*wz^2*l0^2*PotentialPenning(u) + (m*w-0.5*q*B)*sum(vx.*y - vy.*x); % - minEnergy;
% %      PlanarKEnergy_rot(i) = 0.5*m*sum(vx.^2 + vy.^2);
% %      Wall(i) = sum(x.^2-y.^2);
% %      rotational(i) = 0.5*m*w^2*sum(x.^2+y.^2);
% %     %PlanarEnergy(i) =  PlanarKEnergy(i)+PlanarPotEnergy(i);
%    
%     
% end
% 
% plot(AxialEnergy)
% hold on
% plot(PlanarEnergy,'g')
% plot(PlanarEnergy+AxialEnergy,'k')
% plot(AxialKEnergy)
% hold on
 % plot(AxialPotEnergy,'g')
% plot(PlanarKEnergy,'k')
% hold on
% plot(PlanarPotEnergy,'r')

% plot(PlanarKEnergy+PlanarPotEnergy,'k')
% hold on
% plot(PlanarPotEnergy_rot+PlanarKEnergy_rot-rotational,'r')

% semilogy(AxialEnergy)
% hold on
% semilogy(PlanarEnergy,'g')
% semilogy(PlanarEnergy+AxialEnergy,'k')
