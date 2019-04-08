%% Build 2D layer points around a gold nanorod
%  The nanorod is divided into a rod part and two sphere parts. 
%  The rod parts has points with log 
%  Input:  Diameter, height: shape parameters of gold nanorod
%          range: in layer case, the range is the distance of dipole layer
%          points to the surface of particle
%          nrad:  number of angles for the points outside of the sphere
%          parts

function [ x, z ] = RodLayer2D ( diameter, height, range )

%% Build layer range for the rod part
z_mesh = range ; % distance from surface to dipole layers 
x_mesh = round( height/2 - diameter/2 );
x1 = linspace( 0, x_mesh, x_mesh );
for i = 1 : length (x1)
z1(i) = diameter/2 + z_mesh;
end 
%% Inbed points on the sphere part of the nanorods
% The angular interval of points on the sphere part should have the same
% density as that on the side part

theta = linspace(0, pi/2, round( 1.5 * x_mesh/(height/2 - diameter/2) * 0.25 * pi * diameter) );
theta(end) = []; % remove the pi/2 angle because it was covered in the rod parts 
for i = 1 : length(theta)
    [x2(i, :), ~, z2(i, :)] = sph2cart( 0, theta(i), diameter/2 + range); % radial parts has the same number of elements as the range
end
% shift coordinates to the tip of nanorods
x2 = x2 +  height/2 - diameter/2 ;

%% Combine points on rod and sphere parts
x = [x1(:);x2(:)];
z = [z1(:);z2(:)];

end 