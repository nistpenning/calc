function v = spin_down(u,v,f)
    global l0
    N = length(u)/2;
    radii = l0*sqrt(u(1:end/2).^2 + u(end/2+1:end).^2);
    velocities = 2*pi*f*1e3*radii;
    for i = 1:N
        rot = [-u(i+N) u(i)];
        rot = rot./norm(rot);
        % counter velocity to move to rotating frame
        rot = -velocities(i)*rot; 
        %v(i) = rot(1);      
       % v(i+N) = rot(2);

       v(i) = v(i) + rot(1);      
       v(i+N) = v(i+N) + rot(2);
    end
end