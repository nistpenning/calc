% norm_coords = modeBasis(z,E,modes)
%
% Write z positions in normal mode basis <E> for modes <modes>
% 
% Must have trap parameters set by setTrapParameters(...)

function norm_coords = modeBasis(z,E,modes)

N = length(z);

if nargin < 3
    modes = 1:N;
end
  
norm_coords = zeros(1,N);

for m = modes
    norm_coords(m) = dot(z',E(:,m));
end

end

