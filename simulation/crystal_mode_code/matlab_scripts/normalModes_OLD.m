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

function [E,D,st,A,T] = normalModes_OLD(u,zflag)

global N wr wc G Vw M 	                % Load relevant parameters 

[r,xsquare,ysquare] = pairDistance(u);  % Calculate pairwise distances 

%------------------------
%      Axial Modes
%------------------------
if zflag 

    A=zeros(N); %N x N tensor for second derivative at equilibrium positions for axial direction    
    for n=1:N
        A(n,[1:n-1 n+1:N])= 0.5*(r(n,[1:n-1 n+1:N])).^-3;   % alpha != beta, abs() returns complex magnitude
        A(n,n) = 1 - 0.5*sum( (r(n,[1:n-1 n+1:N])).^-3 );   % alpha == beta (diagonal)
        A(n,:) = A(n,:)./(sqrt(M(n)*M));                    % Divide matrix by sqrt(m_i*m_j)
    end
  
    [E,D] = eig(A);    % Calculate eigenvectors E, eigenvalues D
	D = diag(D);       % Extract diagonal (eigenvalues on diagonal of D)
    
    [D ind] = sort(D); % Sort by ascending eigenvalues
    E = E(:,ind);      % Match eigenvectors to sorted eigenvalues
    st = isreal(D);    % Stability criteria. If all of the frequencies are real, the crystal is stable as one plane.
    
    % Convert to real eigenfrequencies and eigenvectors (not orthogonal)
    D = sqrt(D/mean(M));
    E = sqrt(mean(M))*(diag(M.^(-.5))*E); % Convert to original eigenvectors
    T = [];
    
%------------------------
%     Planar Modes
%------------------------
else 

    A=zeros(2*N); %2N x 2N tensor for second derivative at equilibrium positions
    
    % Off diagonal (particle and direction), mass independent
    for i=1:N 
        % Derivatives of x forces (Axx) (Top Left) (alpha !=beta)
        A(i, [1:i-1 i+1:N]) = -0.5*(((r(i,[1:i-1 i+1:N])).^2 - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Derivatives of y forces (Ayy) (Bottom Right) (alpha !=beta)
        A(i+N, [N+1:N+i-1 i+1+N:2*N]) = -0.5*(((r(i,[1:i-1 i+1:N])).^2 - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Mixed derivatives (Axy) (Top Right)
        A(i, [N+1:N+i-1 N+i+1:2*N]) = 0.5*3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5; % alpha != beta
        A(i, N+i) = -0.5*3*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                 % alpha == beta
        
        % Mixed derivatives (Ayx) (Bottom Left)
        A(i+N, [1:i-1 i+1:N]) = 0.5*3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5;       % alpha != beta
        A(i+N, i) = -0.5*3*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                 % alpha == beta           
       
    end
  
    for i=1:N
        % Rewrite diagonal (alpha == beta)
        A(i,i) = (wr^2 - wr*wc + 0.5 - G*Vw) ... 
            + 0.5*sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        A(i+N,i+N) = (wr^2 - wr*wc + 0.5 + G*Vw)  ...
            + 0.5*sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
    end
    
    %Identity Matrix (2N x 2N)
    I = eye(2*N); 
    
    %Zero Matrices 
    Z_2N = zeros(2*N); %(2N x 2N)
    Z_N  = zeros(N); %(N x N)

    % T Matrix 
    T =  wc*([[Z_N eye(N)]; [-eye(N) Z_N]]);
    T = [[Z_N 2*wr*eye(N)]; [-2*wr*eye(N) Z_N]] - T;
    
   % size(A)
   % size(T)
   % size(diag(M))
   
    %mp.Digits(34);
    %MannyLinearization = mp([zeros(2*N) eye(2*N); -V/2 A]); %multiprecision toolbox!
    
   
    %Matlab's polyeig (general eigenvalue problem)
    [E,D] = polyeig(-A,-1i*T,-diag([M M]));
   % size(E)
    D = real(D);
    gz = (D>0);
    E = E(:,gz);
    D = D(gz);
    [D ind] = sort(D);
    %Eraw = E(:,ind); %sorted eigenvectors
    E = E(:,ind); %sorted eigenvectors
    %E=abs(E);
    %E=real(Eraw);
    st = (length(D) == 2*N);


    
    %     for i=1:N
%        A(i,:) =  A(i,:)./M(i);      % Not sure how to deal with defects for planar modes yet
%        A(i+N,:) =  A(i+N,:)./M(i);
%     end
 
%         % T Matrix 
%     T =  wc*([[Z_N eye(N)]; [-eye(N) Z_N]]);
%     T = [[Z_N 2*wr*eye(N)]; [-2*wr*eye(N) Z_N]] - T;
%     
%     %Matlab's polyeig (general eigenvalue problem)
%     [E,D] = polyeig(-A,-1i*T,-I);
%     D = real(D);
%     gz = (D>0);
%     E = E(:,gz);
%     D = D(gz);
%     [D ind] = sort(D);
%     %Eraw = E(:,ind); %sorted eigenvectors
%     E = E(:,ind); %sorted eigenvectors
%     %E=abs(E);
%     %E=real(Eraw);
%     st = (length(D) == 2*N);
end

end

