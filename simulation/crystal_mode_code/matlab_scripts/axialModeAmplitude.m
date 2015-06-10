% axialModeAmplitude(FileLocation,modes,binsize)
%
% Project each time step into axial normal coordinates
% and fit time slices to sine wave to extract amplitude to find
% energy of mode
%
% saves amps and freqs to file

function axialModeAmplitude(FileLocation,modes,binsize)

    thetas = dlmread([FileLocation 'thetas.dat']);
    params = dlmread([FileLocation 'params.dat']);
    setTrapParameters(0,0,0);
    global G N
    setTrapParameters(params(2),-params(3)/G,params(1));

    for k = 1:params(5)/binsize
   
        norm_coords = zeros(params(1),binsize);

        % Use current (every <binsize> steps) configuration to find eigenvectors
        filename = [FileLocation int2str((k-1)*params(5)/binsize) '.dat']; 
        M = dlmread(filename);
        u = convertPythonDataToMatlab(M);
        u = rotate(u,-thetas((k-1)*params(5)/binsize); 
        [E,D,st] = normalModes(u,1);
        
        % Save these frequencies
        dlmwrite([FileLocation 'freqs' num2str(params(5)/binsize) '.dat'],D','-append','precision',10);

        for i = 1:binsize

            filename =[FileLocation int2str((k-1)*params(5)/binsize + i-1) '.dat'];
            M = dlmread(filename);
            z = M(3,:);    
            norm_coords(:,i) = modeBasis(z,E,127)';
            
        end

        amps = zeros(1,N);
        maxamps = zeros(1,N);
        t = linspace(0,binsize*params(6)*params(7),binsize);     % time points
    
 %        plot(t',norm_coords(127,:)')
%         pause(1)
         
        for j = modes
            f = fit(t',norm_coords(j,:)','sin1'); % fit to sine wave
            maxamps(j) = max(abs(norm_coords(j,:)));
            amps(j) = f.a1;
            %clf
%             if j == 127
%             plot(t',norm_coords(j,:)')
%             hold on
%             plot(t',f.a1*sin(f.b1*t'+f.c1),'g')
%             pause(1)
%             end
        end

        % Save these amplitudes
        dlmwrite([FileLocation 'amps' num2str(params(5)/binsize) '.dat'],amps,'-append','precision',10);
        dlmwrite([FileLocation 'maxamps' num2str(params(5)/binsize) '.dat'],maxamps,'-append','precision',10);

    end
end





