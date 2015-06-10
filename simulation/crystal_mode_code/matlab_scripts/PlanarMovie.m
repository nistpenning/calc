% Make Planar Eigenvector Video
cleanupObj = onCleanup(@CLEAN);
setTrapParameters(f,Vw,N); % Initial set to get variables
global N wr wc m q ke B wz V0 w G Vw
load('PlaneTransition','plane_tran');

scrsz = get(0,'ScreenSize');

numspec1=0;
numspec2=0;
masspec1=0;
masspec2=0;

N = 217;
Vw =  (0.04^2)/G;    %\omega_Vw / \omega_z
fmax = plane_tran((plane_tran(:,1) == N),2) - 200/1e3;
wcyc = wc*wz;    %compute cyclotron frequency in real units
fmin0 = (wcyc/2 - sqrt((wcyc^2)/4 - (wz^2)/2))/2/pi/1e3;

fmin = (wcyc/2 - sqrt((wcyc^2)/4 - (wz^2)/2 - 2*q*Vw*G*V0/m))/2/pi;
fmin = (fmin + 200)/1e3;
figure;

set(gcf,'color','w')
%vidpath = 'test_video';
%vwr = VideoWriter(vidpath, 'Motion JPEG AVI');


%open(vwr);

stretch = 30;
mode = 10;
disp('Calculating Normal Modes...')
[junk,junk2,junk3,E_raw] = normalModes(u,0);
mg = abs(E_raw);
ph = angle(E_raw);

for wt = 0 : 2*pi/50 : 10*pi
    iter = 1;
   
    quiver(u(1:N)',u(N+1:2*N)',stretch*mg(1:N,mode).*cos(ph(1:N,mode)+wt),...
        stretch*mg(N+1:2*N,mode).*cos(ph(N+1:2*N,mode)+wt),0,'k','LineWidth',1)
    
    hold on
    scatter(0,0,'r+');
    hold off
    axis equal
    axis([-1.2*max(u),1.2*max(u),-1.2*max(u),1.2*max(u)])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    set(gca,'box','off')
    
    axis off;
    iter = iter+1;

    pause(.00000001)
    %writeVideo(vwr,getframe(gcf));
    %clf
end

%close(vwr);


% us = cell(1,3);
% mg = cell(1,3);
% ph = cell(1,3);
% fcount = 1;
% for f = [fmin,(fmax + fmin0+200/1e3)/2,fmax]
%     setTrapParameters([f,Vw,N]);   
%     foldername = ['Data/' num2str(N) '_' num2str(Vw) '_' num2str(f) '_' num2str(numspec1) '_' num2str(masspec1) '_' ...
%         num2str(numspec2) '_' num2str(masspec2)];
%     
%     load([foldername '/data.mat'], 'u','El')
%     
%     [junk,junk2,junk3,E_raw] = calculateNormalModes(u,0);
%     us{fcount} = u;
%     mg{fcount} = abs(E_raw);
%     ph{fcount} = angle(E_raw);
%     fcount = fcount + 1;
% end
% save(['PAPERFIGURES7_10/planarmoviedata.mat'], 'us','mg','ph')
% load(['PAPERFIGURES7_10/planarmoviedata.mat'], 'us','mg','ph')






