% Look crystal from side (in plane) to see plane transition

% OLD % Range over snapshots to see average z positions of ions

cleanupObj = onCleanup(@CLEAN);

setTrapParameters(45,-80,127);
global l0
FileLocation = 'C:\Users\ACKWinDesk\Desktop\PenningSimulationData\2014_1_27_PlaneTransition\';
vidpath = '1-2PlaneTransition127';
vwr = VideoWriter(vidpath, 'Motion JPEG AVI');
%vwr.FrameRate = 30;
open(vwr);
%thetas = dlmread([FileLocation 'thetas.dat']);
%omegas = dlmread([FileLocation 'omegas.dat']);

for i = 36700:1:37000 %centered around "1-2 PlaneTransition"
    
    for j = 1:10
        %setTrapParameters(omegas(i+1),-80,127);
        filename =[FileLocation int2str(i*10+j) '.dat'];
        M = dlmread(filename);
        [u z] = convertPythonDataToMatlab(M);

        %plot(sort(z))
        %hist(z,100)
        scatter(u(1:end/2),z)
        %scatter(u(end/2+1:end),z)
        axis([-15 15 -.08 .08])

        hold on
        %pause(.001)
    end
    %title([int2str(i) ' ' num2str((46.5e3+i*.02)/1000)])
    title(['Rotation Frequency = ' num2str((46.5e3+i*.02)/1000) ' kHz'])
    xlabel(['Scaled x position, scaling constant l_0 = ' num2str(l0) ' meters'])
    ylabel(['Scaled z position, scaling constant l_0 = ' num2str(l0) ' meters'])
    writeVideo(vwr,getframe(gcf));
    pause(.00001)
    hold off
    
end

close(vwr);

