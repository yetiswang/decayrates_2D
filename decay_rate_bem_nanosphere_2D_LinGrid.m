%% BEM calculation of nanorods

clear all
close all
clc

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  table of dielectric functions
epstab = { epsconst( 1.3^2 ), epstable( 'gold.dat' ) };
%  diameter of sphere
diameter = 60;
Q_int = 0.01 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 41 ),  ...
     pi * linspace( 0, 1, 41 ) .^ 2, diameter, 'triangles' );

% initialize particle
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
    
%%  dipole oscillator
enei = 550;

% 2D positions of dipole
x_mesh = round ( diameter/2 + 50 );
z_mesh = round ( diameter/2 + 50 );

[ x, z ] = meshgrid( linspace( - x_mesh, x_mesh, x_mesh ), linspace( - z_mesh, z_mesh, z_mesh ) );
pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 0.8 );

%  dipole excitation , x , y and z direction
dip = dipole( pt, [ 1, 0, 0 ; 0 , 1 , 0 ; 0, 0, 1 ], op );

figure
bem_plot(p)
hold on
plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3) ,'r.')
hold off
saveas(gcf, 'dipole-sphere.fig')
saveas(gcf, 'dipole-sphere.png')
%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei );
%  total and radiative decay rate
[ tot, rad, rad0 ] = dip.decayrate( sig );

%%  final plot
%% plot 2d map of quantum yield enhancement
QY_x = rad(:,1)./(tot(:,1) + (1- Q_int)/Q_int ) ;
QY_y = rad(:,2)./(tot(:,2) + (1- Q_int)/Q_int ) ;
QY_z = rad(:,3)./(tot(:,3) + (1- Q_int)/Q_int ) ;
QY_average = (rad(:,1) + rad(:, 2) + rad(:, 3))/3 ./( (tot(:,1) + tot(:, 2) + tot(: ,3) )./3 + (1 - Q_int) / Q_int  );

figure
subplot(2,2,1)
scatter3( pt.pos(:,1), pt.pos(:,3), QY_x, 8, QY_x, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('x dipole orientation')
subplot(2,2,2)
scatter3( pt.pos(:,1), pt.pos(:,3), QY_y, 8, QY_y, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('y dipole orientation')
subplot(2,2,3)
scatter3( pt.pos(:,1), pt.pos(:,3), QY_z, 8, QY_z, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('z dipole orientation')
subplot(2,2,4)
scatter3( pt.pos(:,1), pt.pos(:,3), QY_average, 8, QY_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('average dipole orientation')
maximize(gcf)
saveas(gcf, 'QY enhancement_2D.fig')
saveas(gcf, 'QY enhancement_2D.png')
%% plot 2d map decay rates
figure
subplot(2,3,1)
scatter3( pt.pos(:,1), pt.pos(:,3), rad(: , 1 ),8, rad(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('x dipole radiative decay rate')
subplot(2,3,2)
scatter3( pt.pos(:,1), pt.pos(:,3), rad(: , 2 ),8, rad(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('y dipole radiative decay rate')
subplot(2,3,3)
scatter3( pt.pos(:,1), pt.pos(:,3), rad(: , 3 ),8, rad(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('z dipole radiative decay rate')
subplot(2,3,4)
scatter3( pt.pos(:,1), pt.pos(:,3), tot(: , 1 ),8, tot(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('x dipole total decay rate')
subplot(2,3,5)
scatter3( pt.pos(:,1), pt.pos(:,3), tot(: , 2 ),8, tot(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('y dipole total decay rate')
subplot(2,3,6)
scatter3( pt.pos(:,1), pt.pos(:,3), tot(: , 3 ),8, tot(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title('z dipole total decay rate')
maximize(gcf)
saveas(gcf, 'decayrates_2D.fig')
saveas(gcf, 'decayrates_2D.png')
%% save data
dipole_pos = pt.pos ; 
save gamma_rad_BEM_2D rad
save gamma_tot_BEM_2D tot
save dipole_pos_2D dipole_pos
save QY_average_BEM_2D QY_average
save QY_x_2D QY_x
save QY_z_2D QY_z
save QY_y_2D QY_y



