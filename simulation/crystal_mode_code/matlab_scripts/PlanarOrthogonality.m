%addpath('C:\Users\ACKWinDesk\Documents\Multiprecision Computing Toolbox\')
%N = 19;  % Number of ions
%setTrapParameters(54.5,-400,19);
setTrapParameters(64,-80,7);% Set global variables
%N = 127;  % Number of ions
%setTrapParameters(44,-80,N);% Set global variables
global M 
Mass = diag([M M]);
u0 = generateLattice(N,1);  % make triangular lattice with hexagonal boundary
u0 = rotate(u0,pi/2);       % get same major axis as dominics code
u0 = findEquilibrium(u0);    % Solve for equilibrium

% Polyeig (Matlab's polynomail eigenvalue problem solver)
%[Ep,Dp,st,V,A] = normalModes_5_27_14(u,0); % Planar Modes with correction
[Ep,Dp,st,V,A] = normalModes(u0,-1); % Planar Modes with correction
%[Ep,Dp,st,V,A] = normalModes_OLD(u0,0); % Planar Modes with correction

Ep(2*N+1:end,:) = 1i*Ep(1:2*N,:).*repmat(Dp,1,2*N)'; % enforce eigenvectors to obey time derivative

% Linearizations
%MannyLinearization = mp([zeros(2*N) eye(2*N); V/2 A]);
%mp.Digits(34);
MannyLinearization = [zeros(2*N) eye(2*N); V/2 A];
% OtherLinearization = [-1i*A/2 -0.25*A*A + V/2; eye(2*N) -1i*A/2];

% Find eigenvectors and eigenvalues
[EM, DM] = eig(MannyLinearization);
%[EM, DM] = polyeig(V/2,-1i*A,diag([M M]));
% [EO, DO] = eig(OtherLinearization);

% Manny's Eigenvectors
% DM = diag(DM);
% DM = imag(DM); % eigenvalues are now real frequencies
% gz = (DM>0);   % get positive values
% EM = EM(:,gz); % only use those eigenvectors
% DM = DM(gz);
% [DM ind] = sort(DM); % sort in ascending order
% EM = EM(:,ind); % sorted eigenvectors
% 
% % Other Eigenvectors
% DO = diag(DO);
% DO = real(DO); % eigenvalues were already real
% gz = (DO>0);   % get positive values
% EO = EO(:,gz); % only use those eigenvectors
% DO = DO(gz);
% [DO ind] = sort(DO); % sort in ascending order
% EO = EO(:,ind); % sorted eigenvectors

% for i = 1:2*N   % I was wrong, this other linearization is not conventionally orthonormal
%     for j = 1:2*N
%         OtherOrthogonality(i,j) = dot(EO(:,i),EO(:,j));
%     end
% end
% max(max(abs(OtherOrthogonality-(diag(ones(1,2*N)))))) % shows for largest absolute overlap (off diagonal terms only)
% 
% 
for i = 1:2*N   % I was wrong, this other linearization is not conventionally orthonormal
    for j = 1:2*N
       %EnergyM(i,j) = imag(EM(1:2*N,i))'*(Mass*(DM(i)*DM(j)))*imag(EM(1:2*N,j)) - 0.5*real(EM(1:2*N,i))'*V*real(EM(1:2*N,j));
       %EnergyM(i,j) = (EM(1:2*N,i))'*((DM(i)*DM(j)) - 0.5*V)*(EM(1:2*N,j));   %THIS DOESN'T WORK?
       %EnergyM(i,j) = (EM(1:2*N,i))'*(DM(i)*DM(j))*(EM(1:2*N,j))  -0.5*EM(1:2*N,i)'*V*EM(1:2*N,j);
       %EnergyM(i,j) = (EM(1:2*N,i))'*(DM(i)*DM(j))*(EM(1:2*N,j))  - EM(1:2*N,i)'*V*EM(1:2*N,j);
       % KEM(i,j) = (EM(1:2*N,i))'*(DM(i)*DM(j))*(EM(1:2*N,j));
       % PEM(i,j) = -0.5*EM(1:2*N,i)'*V*EM(1:2*N,j);
       EnergyP(i,j) = (Ep(1:2*N,i))'*(Dp(i)*Dp(j))*(Ep(1:2*N,j))  -0.5*Ep(1:2*N,i)'*V*Ep(1:2*N,j);
       %EnergyP(i,j) = 0.5*(Ep(1:2*N,i))'*(Dp(i)*Dp(j))*(Ep(1:2*N,j))  - 0.5*Ep(1:2*N,i)'*V*Ep(1:2*N,j);
       %KEP(i,j) = (Ep(1:2*N,i))'*(Dp(i)*Dp(j))*(Ep(1:2*N,j));
       % PEP(i,j) = -0.5*Ep(1:2*N,i)'*V*Ep(1:2*N,j);
       %EnergyM(i,j) = (EM(2*N+1:end,i))'*(EM(2*N+1:end,j)) - 0.5*EM(1:2*N,i)'*V*EM(1:2*N,j);
       
       %EnergyO(i,j) = (EO(2*N+1:end,i))'*(EO(2*N+1:end,j)) - EO(1:2*N,i)'*V*EO(1:2*N,j);
       %EnergyO(i,j) = (EO(1:2*N,i))'*((DM(i)*DM(j)) -0.5*V)*(EO(1:2*N,j));
       
       %Energy(i,j) = imag(EO(1:2*N,i))'*(Mass/(DO(i)*DO(j)))*imag(EO(1:2*N,j)) - 0.5*real(EO(1:2*N,i))'*V*real(EO(1:2*N,j));
    end
end
% 
% for i = 1:2*N 
%     DMS(i) = mean(imag(EM(2*N+1:end,i)./EM(1:2*N,i)));
% end
% 
% 
% Mass = diag([M M]);
% disp('------------')
% mode2=1;
% mode1=1;
% Energyp = imag(Ep(:,mode1))'*(Mass/(Dp(mode1)*Dp(mode2)))*imag(Ep(:,mode2)) - 0.5*real(Ep(1:2*N,mode1))'*A*real(Ep(1:2*N,mode2));
% for mode1 = 1:2*N
%     %mode2=mode1;
% Energy1 = E1(1:2*N,mode1)'*(Mass/D1(mode1)^2 - A)*E1(1:2*N,mode1)/4;
% 
% 
% mg = abs(E1);
% ph = angle(E1);
% scatter(u(1:N)'+stretch*mg(1:N,mode1).*cos(ph(1:N,mode1)+wt),u(N+1:2*N)'+stretch*mg(N+1:2*N,mode1).*cos(ph(N+1:2*N,mode1)+wt),'go', 'filled');