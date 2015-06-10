
loadflag = 1;
if loadflag
    FileLocation = 'D:\PenningSimulationData\2014_10_28_PlanarTemp\';  
    setTrapParameters(0,0,0);
    global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 params
    params = dlmread([FileLocation 'params.dat']);
    setTrapParameters(params(2),-params(3)/G,params(1));
    thetas = dlmread([FileLocation 'thetas.dat']);
    [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,-1);
    for i = 0:params(5)-1
        us(i+1,:) = rotate(us(i+1,:),-thetas(i+1));           
    end
load([FileLocation 'DataVars.mat'],'PKE','PPE')
load([FileLocation 'planarModeDecomposition.mat'],'norm_coords_planar','EnergyMode')
totalE = PPE+PKE';
end

u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium


%us = us(end/2:end,:);
[avg, R, avgdist] = averageDistance(u0); % get lattice scaling with avgdist

[uadd x y] = addHexShell(3);
uadd = [uadd(1:2:end) uadd(2:2:end)];
uadd = rotate(uadd,pi/2);       % get same major axis as dominics code
uadd = avgdist*uadd;
 
utest = [u0(1:N) uadd(1:end/2) u0(1+N:2*N) uadd(end/2+1:end)];

[verts cells] = voronoin([utest(1:end/2)' utest(end/2+1:end)']);
[ushift sigma] = normal_traj(us,u0);

figure
cmp = colormap(hsv(N));

voronoi(u0(1:end/2),u0(end/2+1:end),'k');
%voronoi(utest(1:end/2),utest(end/2+1:end));
hold on
%scatter(ushift(1:end/2),ushift(end/2+1:end),'filled','w');
scatter(u0(1:end/2),u0(end/2+1:end),'filled','k');

experiments = 12;
cmp2 = colormap(jet(experiments));
for k = 1:experiments
    ind = (4e4*(2*k)-4e4):4e4*(2*k);
    mtE(k) = mean(totalE(ind));
end
[junk, index] = sort(mtE);

for i=1:N
    
    %h = patch(verts(cells{i}',1),verts(cells{i}',2),cmp(i,:));
    h = patch(verts(cells{i}',1),verts(cells{i}',2),'w');
    %[elipx,elipy] = ellipse(2*sigma(i),2*sigma(i+N),u0(i),u0(i+N),cmp(i,:),'k',1);
    
    for k = 1:experiments % cooling lasers off
        ind = (4e4*(2*index(k))-4e4):4e4*(2*index(k));
        [umean sigma] = normal_traj(us(ind,:),u0);
        ellipse(2*sigma(i),2*sigma(i+N),u0(i),u0(i+N),'w',cmp2(index(k),:),2);
    end
    
    %ind = ~inpolygon(us(:,i),us(:,i+N),verts(cells{i}',1),verts(cells{i}',2));
    %ind = ~inpolygon(us(:,i),us(:,i+N),elipx,elipy);
    %temp = ~inpolygon(us(40000:80000,i),us(40000:80000,i+N),elipx,elipy);
    %ind = zeros(length(us),1);
    %ind(40000:80000) = temp;
    %plot(us(find(ind),i),us(find(ind),i+N),'k.');
    %plot(us(1:1000:end,i),us(1:1000:end,i+N),'k.');
    %plot(us(40000:80000,i),us(40000:80000,i+N),'w.');
    if all(inpolygon(us(:,i),us(:,i+N),verts(cells{i}',1),verts(cells{i}',2)));
       %ellipse(2*sigma(i),2*sigma(i+N),u0(i),u0(i+N),cmp(i,:),'w',1)
    else 
       %ellipse(2*sigma(i),2*sigma(i+N),u0(i),u0(i+N),cmp(i,:),'k',1)
        scatter(u0(i),u0(i+N),'filled','m');
    end
end



  
  
  