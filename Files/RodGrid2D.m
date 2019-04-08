%% Build 2D grid points around a gold nanorod
%  The nanorod is divided into a rod part and two sphere parts. 
%  The rod parts has points with log 
%  input: diameter and height of GNR
%         range: spherical region of interst around the GNR
%         nrad:  number of points on the radial part
%         scale: by default the lateral grid is per nm, use scale to
%         increase of decrease the density.
function [ x, z ] = RodGrid2D ( diameter, height, range, nrad, scale )

if nargin == 4 
    scale = 1;
else 
end 
%% Build grid range for the rod part
z_mesh = round( ( diameter/2 + range));
x_mesh = round( (height/2 - diameter/2) );
[ x1, z1 ] = meshgrid( linspace( 0 , x_mesh,  x_mesh.*scale  ), logspace( log10(z_mesh - range), log10(z_mesh),  nrad ) );

%% Inbed points on the sphere part of the nanorods
% The angular interval of points on the sphere part should have the same
% density as that on the side part
theta = linspace(0, pi/2, round( 1.5 * x_mesh/(height/2 - diameter/2) * 0.25 * pi * diameter) .*scale );

theta(end) = []; % remove the pi/2 angle because it was covered in the rod parts 
for i = 1 : length(theta)
    [x2(i, :), ~, z2(i, :)] = sph2cart( 0, theta(i), logspace( log10(diameter/2 ) ,log10(diameter/2 + range), nrad )); % radial parts has the same number of elements as the range
end
% shift coordinates to the tip of nanorods
x2 = x2 +  height/2 - diameter/2 ;

%% Combine points on rod and sphere parts
x = [x1(:);x2(:)];
z = [z1(:);z2(:)];

end 