
function [ushift sigma] = normal_traj(us_rot,u0)

    motion = us_rot - repmat(u0,length(us_rot),1); % subtract off equilibrium positions
    [mu,sigma] = normfit(motion);
    ushift = mu + u0;

end