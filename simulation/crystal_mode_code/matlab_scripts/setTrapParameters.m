% setTrapParameters(f,vw,n)
%
% Declares trap parameters as global, generates mass array
%
% INPUT 
% f:  Crystal/Rotating Wall rotation frequency in kHz,
% vw: Potential applied to quadrupole "rotating wall" (peak-to-peak/2)
% n:  Number of ions (N)
%
% EXAMPLE: 
% setTrapParameters(43,105,217) 
% 	43 kHz rotation frequency
% 	105 V for "strong rotating wall"
% 	N = 217 ions

function setTrapParameters(f,vw,n)

%Declare as global variables, every function that also needs these
%parameters must have this same line
global N m wr wc q ke B wz V0 w G Vw ww M l0 v0 t0

% Fixed constants (MKS units)
mu = 1.660538921e-27;         % mass unit
m  = 9.012182*mu;             % atomic mass of berylium (ion negligible difference)
q  = 1.602176565*10^(-19);    % electron charge
ke = 8.987551787368*10^9;     % coulomb constant
B  = 4.4584;                  % Magnetic Field magnitude (tesla)
V_T = 1000;                   % Potential on endcap 
G = .045/V_T;                 % "Geometric factor" for "rotating wall" (relates V0 to Vwall)
wz = 2*pi*795e3;              % Axial Frequency of Penning trap
V0 = (0.5*m*wz^2)/q;          % Effective Axial potential

% Tunable Parameters
w  = 2*pi*f*1e3;              % Quadrupole "rotating wall" frequency (same rotation frequency as crystal)
Vw = vw;                      % Actual potential on rotating wall 
ww = sqrt(2*q*Vw/m);          % Effective trapping frequency of rotating wall
N = n;                        % Number of ions

% Convert to dimensionless quantities
l0 = ((ke*q^2)/(q*V0))^(1/3); % Dimensionless length 
t0 = 1/wz;                    % Dimensionless time
v0 = l0/t0;                   % Dimensionless velocity
wr = w/wz;       			  % Dimensionless crystal rotation frequency
wc = q*B/(m*wz); 			  % Dimensionless cyclotron frequency

% Dimensionless mass array
M = ones(1,N);                % All ions are same species where m is scaled to 1
							  % To add defects see: addDefects.m
end



