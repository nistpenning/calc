% u = addDefects(u0,numspec1,masspec1,numspec2,masspec2,numspec0,masspec0)
%
% Iteratively replaces original ions with defects. Here, original == Be^+
%
% First, replaces heavy defects (species 2 then species 1)
% At each iteration an original ion with the largest radius from the origin 
% is replaced and the the equilbrium structure is re-solved.
% 
% Second, replaces light defects (species 0)
% At each iteration an originan ion with the smallest radius from the origin 
% is replaced and the the equilbrium structure is r-esolved.
%
% This method is motivated by the effect of the centrigual force on defects
% There is probably a more elegant way to do this for an arbitrary number
% of defect types but this works fine. It's slow because each iteration
% it recalculates the equilbrium position
%
% If no ions need to be replaced, the functions skips to the last line
% and computes the equilbrium lattice for a homogeneous crystal
%
% INPUT 
% u0:        Seed lattice
% numspec1:  number of species 1 (heavy) ions (still have elementary charge)
% masspec1:  mass of species 1 (heavy) ions relative to original mass (Dimensionless mass scaled by original mass)
% numspec2:  number of species 2 (heavy) ions (still have elementary charge)
% masspec2:  mass of species 2 (heavy) ions relative to original mass (Dimensionless mass scaled by original mass)
% numspec0:  number of species 0 (light) ions (still have elementary charge)
% masspec0:  mass of species 0 (light) ions relative to original mass (Dimensionless mass scaled by original mass)
%
% OUTPUT
% u:         Equilibrium lattice with defects
%
% NOTE:      "Identity" of ion is stored in mass array M ("identity" is dimensionless mass where 1 is original mass)
%
% EXAMPLE
% For N = 217, u0 is output from generateLattice.m or idealLatticeScale.m
% u = addDefects(u0,24,1.128,36,2.8863,1,0.1118);
%
% This adds 36 BeOH ions with scaled mass 2.8863 = (9.012182+17)/9.012182
%      and  24 BeH  ions with scaled mass 1.1280  = (9.012182+1.00794)/9.012182
%      and  1  H    ion  with scaled mass 0.1118 = (1.00794)/9.012182 
%###################################################################################################################

function u = addDefects(u,numspec1,masspec1,numspec2,masspec2,numspec0,masspec0)

global N M   % Load relevant parameters (must have called setTrapParameters(...))

% Checks if defects are to be added
if numspec0 || numspec1 || numspec2 > 0
    
    if numspec0 + numspec1 + numspec2 > N
        disp('Error: More defects than number of ions!')
        return
    end

    % Replace with Species 2
    for outer = 1:numspec2
        R = sqrt(u(1:N).^2 + u(N+1:2*N).^2);   % Find each ion's distance from origin
        Be = find(M==1);                       % Indices that are original mass
        R = R(Be);                             % Just look at radius of those indices find radius 
        M(Be(find(R == max(R),1))) = masspec2; % Set first max radius ion's mass to mass of species 2
        [u, V] = findEquilibrium(u);           % Solve new system with seed lattice as solution from previous iteration
    end
    
    % Replace with Species 1
    for inner = 1:numspec1
        R = sqrt(u(1:N).^2 + u(N+1:2*N).^2);   % Find each ion's distance from origin      
        Be = find(M==1);                       % Indices that are original mass
        R=R(Be);                               % Just look at radius of those indices find radius 
        M(Be(find(R == max(R),1))) = masspec1; % Set first max radius ion's mass to mass of species 1
        [u, V] = findEquilibrium(u);           % Solve new system with seed lattice as solution from previous iteration
    end

    % Replace with Species 0
    for inner = 1:numspec0
        R = sqrt(u(1:N).^2 + u(N+1:2*N).^2);   % Find each ion's distance from origin      
        Be = find(M==1);                       % Indices that are original mass
        R=R(Be);                               % Just look at radius of those indices find radius 
        M(Be(find(R == min(R),1))) = masspec0; % Set first min radius ion's mass to mass of species 0
        [u, V] = findEquilibrium(u);           % Solve new system with seed lattice as solution from previous iteration
    end
    
% No defects to add - solve homogeneous system
else 
    [u, V] = findEquilibrium(u);
end

end