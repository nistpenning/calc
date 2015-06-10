
Vwalls = [0,-35,-80,-105];
for i=5:10
    Ns(i-4) = 1 + 6*sum(1:i); % Number of ions in configuration
end


%plane_tran = zeros(length(Ns),length(Vwalls));
%maxf = 49;

%fsres = 0.1; % search resolution, kHz

load('PlaneTransition_2014_1_24','plane_tran')
fsres = 0.01; % search resolution, kHz

for i = 1:length(Vwalls)
    Vw = Vwalls(i);
    %f_start = maxf;
    for j = 1:length(Ns)
        f_start = plane_tran(i,j)+.11;
        N = Ns(j);

        u0 = generateLattice(N,1);
        f_next = f_start;
        setTrapParameters(f_next,Vw,N);
        u = findEquilibrium(u0);
        [E,D,st] = normalModes(u,1);
        while ~isreal(D)   % search down until find stable config

            f_next = f_next - fsres;
            setTrapParameters(f_next,Vw,N);
            u = findEquilibrium(u0);
            [E,D,st] = normalModes(u,1);

        end       
        
        f_start = f_next;
        plane_tran(i,j) = f_start;
        disp([Vw,N,f_start])
        save('PlaneTransition_2014_1_24_01','plane_tran')
    end
end


