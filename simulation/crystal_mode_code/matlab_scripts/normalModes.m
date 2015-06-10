% [E,D,st,r] = normalModes(u,zflag)
%
% Calculates axial or planar eigenvectors and eigenvalues
% Planar eigenvectors with defects has NOT been worked out yet
%
% INPUT
% u:     Equilibrium lattice
% zflag: 1 for axial modes, 0 planar modes
%
% OUTPUT
% E:     Matrix of eigenvectors (columns are eigenvectors)
% D:     List of respective frequencies  
% st:    Stability flag
%		 	1-2 plane transition for axial modes
% 			Planar stability for planar modes 
%			1 == stable, 0 == unstable
%
% Must have trap parameters set by setTrapParameters(...)
%
% Needs to be slightly modified for different masses


function [E,D,st,V,A] = normalModes(u,zflag)

global N wr wc G Vw M 	                % Load relevant parameters 

[r,xsquare,ysquare] = pairDistance(u);  % Calculate pairwise distances 

%------------------------
%      Axial Modes
%------------------------
if zflag == 1

    V=zeros(N); %N x N tensor for second derivative at equilibrium positions for axial direction    
    for n=1:N
        V(n,[1:n-1 n+1:N])= 0.5*(r(n,[1:n-1 n+1:N])).^-3;   % alpha != beta, abs() returns complex magnitude
        V(n,n) = 1 - 0.5*sum( (r(n,[1:n-1 n+1:N])).^-3 );   % alpha == beta (diagonal)
        V(n,:) = V(n,:)./(sqrt(M(n)*M));                    % Divide matrix by sqrt(m_i*m_j)
    end
  
    [E,D] = eig(V);    % Calculate eigenvectors E, eigenvalues D
	D = diag(D);       % Extract diagonal (eigenvalues on diagonal of D)
    
    [D ind] = sort(D); % Sort by ascending eigenvalues
    E = E(:,ind);      % Match eigenvectors to sorted eigenvalues
    st = isreal(D);    % Stability criteria. If all of the frequencies are real, the crystal is stable as one plane.
    
    % Convert to real eigenfrequencies and eigenvectors (not orthogonal)
    D = sqrt(D/mean(M));
    E = sqrt(mean(M))*(diag(M.^(-.5))*E); % Convert to original eigenvectors
    A = [];
    
%------------------------
%     Planar Modes
%------------------------
else 

    V=zeros(2*N); %2N x 2N tensor for second derivative at equilibrium positions
    
    % Off diagonal (particle and direction), mass independent
    for i=1:N 
        % Derivatives of x forces (Axx) (Top Left) (alpha !=beta)
        V(i, [1:i-1 i+1:N]) = -(((r(i,[1:i-1 i+1:N])).^2 - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Derivatives of y forces (Ayy) (Bottom Right) (alpha !=beta)
        V(i+N, [N+1:N+i-1 i+1+N:2*N]) = -(((r(i,[1:i-1 i+1:N])).^2 - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Mixed derivatives (Axy) (Top Right)
        V(i, [N+1:N+i-1 N+i+1:2*N]) = 3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5; % alpha != beta
        V(i, N+i) = -3*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                 % alpha == beta
        
        % Mixed derivatives (Ayx) (Bottom Left)
        V(i+N, [1:i-1 i+1:N]) = 3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5;       % alpha != beta
        V(i+N, i) = -3*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                 % alpha == beta           
       
    end
  
    for i=1:N
        % Rewrite diagonal (alpha == beta)
        V(i,i) = 2*M(i)*(wr^2 - wr*wc + 0.5 - G*Vw) ... 
            + sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        V(i+N,i+N) = 2*M(i)*(wr^2 - wr*wc + 0.5 + G*Vw)  ...
            + sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
    end
    
    %Identity Matrix (2N x 2N)
    I = eye(2*N); 
    
    %Zero Matrices 
    Z_2N = zeros(2*N); %(2N x 2N)
    Z_N  = zeros(N); %(N x N)

    % T Matrix 
    A =  wc*([[Z_N eye(N)]; [-eye(N) Z_N]]);
    A = [[Z_N 2*wr*eye(N)]; [-2*wr*eye(N) Z_N]] - A;
    
    % Linearizations
    %mp.Digits(34);
    %MannyLinearization = mp([zeros(2*N) eye(2*N); V/2 A]); %multiprecision toolbox apparently only a trial version..
    MannyLinearization = [zeros(2*N) eye(2*N); V/2 A]; 
    dlmwrite('MannyLinearization.txt',MannyLinearization,'delimiter',' ','precision',16);
    %!"C:\Program Files\Wolfram Research\Mathematica\8.0\MathKernel" -noprompt -run "<<C:\Users\ACKWinDesk\Documents\GitHub\ultracold-ions\matlab_scripts\MultiprecisionEigensystem.m"    
    !"C:\Program Files\Wolfram Research\Mathematica\8.0\math" -noprompt -run "<<C:\Users\ACKWinDesk\Documents\GitHub\ultracold-ions\matlab_scripts\MultiprecisionEigensystem.m"    

    MathEigValReal = dlmread('MathematicaEigenvalues_real.txt');
    MathEigValImag = dlmread('MathematicaEigenvalues_imag.txt');
    MathEigVecReal = dlmread('MathematicaEigenvectors_real.txt');
    MathEigVecImag = dlmread('MathematicaEigenvectors_imag.txt');

    D = MathEigValReal + 1i*MathEigValImag;
    E = MathEigVecReal' + 1i*MathEigVecImag';
    st = 1;
    %return
    [junk ind] = sort(abs(D)); % sort in ascending order
    D = imag(D(ind));
    E = E(:,ind); % sorted eigenvectors

    % Normalize by energy, rather than norm = 1
   for i = 1:4*N   
       E(:,i) = E(:,i)/sqrt(real((E(1:2*N,i))'*(D(i)*D(i))*(E(1:2*N,i))  -0.5*E(1:2*N,i)'*V*E(1:2*N,i)));
   end

    
    
    % Use Matlab to find eigenvectors and eigenvalues (precision issues?)
    %[E, D] = eig(MannyLinearization);
    %condeig(MannyLinearization)
    %[E,D] = polyeig(V/2,-1i*A,diag([M M]));
    %size(D)
   
    
     %D = diag(D);
     %D = imag(D); % eigenvalues are now real frequencies
     %[D ind] = sort(abs(D)); % sort in ascending order
     %[junk ind] = sort(abs(D)); % sort in ascending order
     %D = D(ind);
     %E = E(:,ind); % sorted eigenvectors
     
    if zflag == 0
        % Cleanup and Sort Eigenvectors
       
        %D = real(D);
        gz = (D>0);   % get positive values
        E = E(:,gz); % only use those eigenvectors
        D = D(gz);
       
    end
   
    %EnergyM(i,j) = imag(EM(1:2*N,i))'*(Mass*(DM(i)*DM(j)))*imag(EM(1:2*N,j)) - 0.5*real(EM(1:2*N,i))'*V*real(EM(1:2*N,j));
  
end

end

