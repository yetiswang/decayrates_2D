%% Calculate the volume of a GNR

function [AR, vol] = GNRVol(height, diameter)

AR = height/diameter; 

vol = 3/4 * pi * ( 0.5*diameter )^3 +  pi * ( 0.5*diameter )^2 * ( height - diameter ) ;

end 
