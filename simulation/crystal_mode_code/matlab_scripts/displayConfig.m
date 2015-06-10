% displayConfig(u)
%
% displays configuration u, returns handle
% black is original ions, (Be)
% red is origin
% magenta is species 0 ions (H)
% green is species 1 ions   (BeH)
% blue is species 2 ions    (BeOH)
%
% must call setTrapParameters before calling this function

function h = displayConfig(u)

global M N

beh = find(M>1);
if ~isempty(beh)
beoh = find(M>min(M(beh)));
beh = setdiff(beh,beoh);
hyd = find(M<1);
else
    beoh = [];
    beh = [];
    hyd = [];
end

%h = figure;

scatter(u(1:end/2),u(end/2+1:end),'ko', 'filled');
hold on
scatter(0,0,'r+');
scatter(u(beh),u(beh + N),'go', 'filled');
scatter(u(beoh),u(beoh + N),'bo', 'filled');
scatter(u(hyd),u(hyd + N),'mo', 'filled');
box on
%scatter(u(beoh),u(beoh + N),'bo', 'filled');
%axis([-1.2*max(u), 1.2*max(u), -1.2*max(u), 1.2*max(u)])
%xlabel('Characteristic Length Units')
%title('Equilibrium Positions')

%hold off
%axis equal
%axis([-240 240 -240 240])
%axis square
%title('Red = Origin, Black = Be, Green = BeH, Blue = BeOH')

end

