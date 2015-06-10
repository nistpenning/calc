% Plot extracted amplitudes for each mode over the time bin
% to see how amplitude for each mode changes over time

amps = dlmread(['ExtractedAmps_2014_2_3_2.dat']);
freqs = dlmread(['ExtractedAmps_2014_2_3_2.dat']);
    
for i = 1:1000
    for j = 1:127
        %semilogy([j j],[1e-12 abs(amps(i,j))])
        plot([j j],[0 abs(amps(i,j))])
        hold on
    end
    %axis([1 127 1e-12 1e-6])
    axis([1 127 0 5e-6])
    pause(.01)
    hold off
end