% Convert python data file and to format for matlab scripts
%
% Read parameters first and call setTrapParameters()

function [us,zs,vxs,vys,vzs] = convertPythonDataToMatlab2(FileLocation,mode)
global params l0
if mode == -1
    M = dlmread([FileLocation 'Data.dat']);
else
    M = dlmread([FileLocation 'Data' num2str(mode) '.dat']);
end

Data = M';
%size(Data)
Data = reshape(Data,params(1),6,params(5));
us = [Data(:,1,:); Data(:,2,:)];
us = squeeze(us)'/l0;

zs = squeeze(Data(:,3,:))';
vxs = squeeze(Data(:,4,:))';
vys = squeeze(Data(:,5,:))';
vzs = squeeze(Data(:,6,:))';

end

%Data = reshape(M,6,params(1),params(5));
%us = squeeze([Data(1,:,:) Data(2,:,:)])'/l0; % reshape to get rows of positions
%zs = squeeze(Data(3,:,:))';
