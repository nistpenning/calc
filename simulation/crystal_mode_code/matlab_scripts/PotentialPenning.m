% [V, F, H, r] = PotentialPenning(u)
%
% Calculates the potential for the collection of ions with positions u
% Format: [x1 x2 ... xN y1 y2 ... yN]
% 
% INPUT
% u:    lattice
% V:    Total (velocity indepednent) potential energy
% F:    Negative force (postive gradient of V) on each ion in each direction [fx(1) fx(2) ... fx(N) fy(1) fy(2) ... fy(N)]
% H:    Hessian matrix (second derivative of V) (positive gradient of F) 
%
% Use with u = fminunc(@PotentialPenning,u0); (see findEquilbrium(u0))
% where u0 is the seed lattice and u is final equilibrium lattice
%
% Trap parameters must be set by setTrapParameters(...)

function [V, F, H, r] = PotentialPenning(u) %, xsquare, ysquare

global N wr wc Vw M G
  
    F = zeros(2*N,1);  % Force (in each direction and each ion)
    H = zeros(2*N);    % Hessian Matrix
    	
    [r,xsquare,ysquare,x,y] = pairDistance(u); % Calculate pairwise distances
    
	%----------------------------------------------------------------------
    % Find Total Potential Energy 
    %----------------------------------------------------------------------
	
    V=1./r;                                  % Find Coulomb potential between ions
    V(isnan(V) | isinf(V))=0;                % Self potential = 0
    %0.5*sum(sum(V))
    V =  -sum( (M.*wr^2 - wr*wc + 1/2).*(x.^2 + y.^2) )  + G*Vw.*sum( (x.^2 - y.^2)  ) + 0.5*sum(sum(V)); % Total potential energy
    
    %----------------------------------------------------------------------
    % Find Total Force (indluding Trap Potentials) (positive gradient of V)
    %----------------------------------------------------------------------
	
	% Find magnitude of coulomb force
    Fc = r.^-2;
    Fc(isnan(Fc) | isinf(Fc)) = 0;  % Self force = 0
       
    % Convert coloumb total forces into cartesian forces fx, fy
    f_x = (xsquare./r).*Fc;         % x component of force
	f_x(isnan(f_x) | isinf(f_x))=0; % Self force = 0
    f_y = (ysquare./r).*Fc;         % y component of force
	f_y(isnan(f_y) | isinf(f_y))=0; % Self force = 0
	
    % Calculate trap potential terms
    Ftrapx = -2*(M.*wr^2 - wr*wc + 1/2 - G*Vw).*x;
    Ftrapy = -2*(M.*wr^2 - wr*wc + 1/2 + G*Vw).*y;
	
	% Add Coulomb and trap forces together
    F(1:N) = -sum(f_x,2) + Ftrapx';
    F(N+1:2*N) = -sum(f_y,2) + Ftrapy';
        
    %----------------------------------------------------------------------
    % Find Hessian Matrix (or Jacobian Matrix, depending on your language)
    %----------------------------------------------------------------------
    
	% Off diagonal (particle and direction), mass independent
    for i=1:N
        % Derivatives of x forces (Axx) (alpha !=beta)
        H(i, [1:i-1 i+1:N]) = (((r(i,[1:i-1 i+1:N])).^2 ...
            - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Derivatives of y forces (Ayy) (alpha !=beta)
        H(i+N, [N+1:N+i-1 i+1+N:2*N]) = (((r(i,[1:i-1 i+1:N])).^2 ...
            - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        
        % Mixed derivatives (Axy)
        H(i, [N+1:N+i-1 N+i+1:2*N]) = - 3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5; % alpha != beta
        H(i, N+i) = 3*0.5*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                % alpha == beta
        
        % Mixed derivatives (Ayx)
        H(i+N, [1:i-1 i+1:N]) = - 3*(xsquare(i,[1:i-1 i+1:N])).*(ysquare(i,[1:i-1 i+1:N])).*(r(i,[1:i-1 i+1:N])).^-5;       % alpha != beta
        H(i+N, i) = 3*sum(xsquare(i,[1:i-1 i+1:N]).*ysquare(i,[1:i-1 i+1:N]).*(r(i,[1:i-1 i+1:N])).^-5);                    % alpha == beta
    end
     
    for i=1:N
        % Rewrite diagonal (alpha == beta)
        H(i,i) = -2*(M(i)*wr^2 - wr*wc + 1/2 - G*Vw) ...
            - sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(xsquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
        H(i+N,i+N) = -2*(M(i)*wr^2 - wr*wc + 1/2 + G*Vw) ...
            - sum( ((r(i,[1:i-1 i+1:N])).^2 - 3*(ysquare(i,[1:i-1 i+1:N]).^2)).*(r(i,[1:i-1 i+1:N])).^-5);
    end
    
end

