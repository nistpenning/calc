N = 217; % S = 8 complete hex shells
setTrapParameters(0,0,N);
global N m wr wc q ke B wz V0 w G Vw ww M
load('PlaneTransition.mat')
fmax = plane_tran((plane_tran(:,1) == N),2) - 200/1e3;
wcyc = wc*wz;    %computer cyclotron frequency in real units
fmin0 = (wcyc/2 - sqrt((wcyc^2)/4 - (wz^2)/2))/2/pi/1e3;
Vw = ([0.04].^2)./G;    %[....] is ratio of \omega_Vw / \omega_z
fmin = (wcyc/2 - sqrt((wcyc^2)/4 - (wz^2)/2 - 2*q*Vw*G*V0/m))/2/pi;
fmin = (fmin + 200)/1e3;

%for f = [fmin,fmax]
for f = [fmax]
    %for numspec1 = [0,1,36]
    for numspec1 = [50,100]
        setTrapParameters(f,Vw,N);
        u0 = generateLattice(N,1);
        u0 = idealLatticeScale(u0,1,0.01);
        u = addDefects(u0,numspec1,1.128,0,2.8863,0,0.1118);
        [Eaxial,Daxial,st_ax] = normalModes(u,1);
        save(['C:\Users\ACKWinDesk\Dropbox\2013_GeorgetownPapers\Defect\DefectData\' num2str(f) '_' num2str(numspec1) '.mat'], 'N' ,'f','wr' ,'wc', 'Vw','M' ,'u','u0','Eaxial','Daxial')
         
    end
end

    
