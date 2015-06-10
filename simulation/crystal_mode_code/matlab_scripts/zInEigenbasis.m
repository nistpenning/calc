% Project each time step into normal coordinates

FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_10_NormalModeExpansionTestWithLaserCooling\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_31_NonEquilibriumTest\';
%FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_2_7_AxialLaserCooling_COM1e-5\';
thetas = dlmread([FileLocation 'thetas.dat']);
params = dlmread([FileLocation 'params.dat']);
setTrapParameters(0,0,0);
global G wz N
setTrapParameters(params(2),-params(3)/G,params(1));


%mode_z = zeros(params(1),params(5));
%for i = 1:1:params(5)
%modes = 45;
%mode45normcoords = zeros(1,1000);
modes = 1:127;
%modes = 31;
binsize = 1000;
for k = 1:params(5)/binsize
%for k = 500:1000
    k
    norm_coords = zeros(params(1),binsize);
    
    % Use current (every <binsize> steps) configuration to find eigenvectors
    filename = [FileLocation int2str((k-1)*params(5)/binsize) '.dat']; 
    M = dlmread(filename);
    [u z] = convertPythonDataToMatlab(M);
    u = rotate(u,-thetas(params(5)-1)); 
    [E,D,st] = normalModes(u,1);
    dlmwrite([FileLocation 'freqs' num2str(params(5)/binsize) '.dat'],D','-append','precision',10);
    for i = 1:binsize

        %filename =[FileLocation int2str(k) '.dat'];
        filename =[FileLocation int2str((k-1)*params(5)/binsize + i-1) '.dat'];
        M = dlmread(filename);
        z = M(3,:);    
        
        norm_coords(:,i) = modeBasis(z,E)';
        %test = modeBasis(z,E,45);
        %mode45normcoords(k) = test(45) ;
    end
    
    amps = zeros(1,N);
    freqs = zeros(1,N);
    t = linspace(0,binsize*params(6)*params(7),binsize);     % time points
%    opts = optimset('Display','off');
    for j = modes
 %       c0 = [max(abs(norm_coords(j,:)))  wz*D(j) 0];
%        c = lsqcurvefit(F,c0,t,norm_coords(j,:),[],[],opts); 
        %f = fit(t',norm_coords(j,:)','sin1','StartPoint',[max(abs(norm_coords(j,:))), wz*D(j),0]);
        %f = fit(t',norm_coords(j,:)','sin1','StartPoint',[max(abs(norm_coords(j,:))), wz*D(j),0], ...
        %    'Lower',[max(abs(norm_coords(j,:)))/2 .9*wz*D(j) -2*pi],'Upper',[1.5*max(abs(norm_coords(j,:))) 1.1*wz*D(j) 2*pi]); % this isn't working right 
        f = fit(t',norm_coords(j,:)','sin1');
        amps(j) = f.a1;
        %freqs(j) = f.b1;
        %phases(j) = f.c1;
    end
   % plot(freqs)
   % pause(1)
%     plot(t,norm_coords(j,:))
%     hold on
%     plot(t,amps(j)*sin(freqs(j)*t+phases(j)),'g')
%     title([num2str(amps(j)) ' ' num2str(freqs(j)) ' ' num2str(phases(j)) ' ' num2str(wz*D(j))])
%     pause(1)
%     hold off

    dlmwrite([FileLocation 'amps' num2str(params(5)/binsize) '.dat'],amps,'-append','precision',10);
   
    
end





