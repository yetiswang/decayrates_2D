%% Build grid points for nanosphere
% inbed points on the sphere part of the nanorods
% The angular interval of points on the sphere part should have the same
% density as that on the side part

function [ x , z ] = SphereGrid2D ( diameter, nangles, range)

theta = linspace( 0, 0.5 * pi, nangles );
for i = 1 : length(theta)
[x(i, :), ~ , z(i,:)] = sph2cart( 0, theta(i), logspace( log10(diameter/2 ) ,log10(diameter/2 + 0.8 + range),range ));
end 

end 