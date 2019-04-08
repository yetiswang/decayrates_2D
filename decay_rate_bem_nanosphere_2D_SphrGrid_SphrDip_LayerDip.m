%% BEM dipole decay rate and quantum yield calculation in the vinicity of nanorods

clear all
close all
clc

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  table of dielectric functions
epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };

%  nanorod geometries
mesh = [ 41, 41, 41]; % n1 for the circumference of the rod, n2 for the polar angles of the rod caps, n3 for the cylinder-shaped middle part of the rod
height =  120   ;
diameter =  38  ;

Q_int = 0.9 ; % intrinsic quantum efficiency
%  nanorod with finer discretization at the top
%  To calculate decay rates close to surface, the mesh close to the
%  positions of dipoles need to refined.
p = trirod ( diameter, height, mesh, 'triangles' );

%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);

%%  dipole oscillator
enei = 800;

% 2D positions of dipole
% Consuming large amount of memory, use mesh size in the units of nanometer 
% inbed points on the side of the nanorods
[ x, z ] = RodLayer2D ( diameter, height, 10);
pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 0.8 );

% calculate the spherical dipole orientation for each dipole location
dir = GetDipDirGNR( height, diameter, pt );

dip = dipoleret( pt, dir, 'full');

%% 
figure
coneplot(pt.pos,dir)
axis on
grid on
%%
saveas(gcf,'dipole-orientation.fig')
saveas(gcf,'dipole-orientation.png')
%%
figure
hold on
bem_plot(p)
plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3) ,'r.')
hold off
view([0 0])
%%
saveas(gcf,'dipole-sphere.fig')
saveas(gcf,'dipole-sphere.png')
%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei );
%  total and radiative decay rate
[ tot, rad, rad0 ] = dip.decayrate( sig );
%% Expand to full angle sphere slice
% Concatenate position of dipoles for the whole 2pi angles
dipole_pos = pt.pos ; 
pos(:,1) = [ dipole_pos(:,1) ; ( - dipole_pos(:,1) ) ; ( - dipole_pos(:,1) ) ; dipole_pos(:,1)];
pos(:,2) = [dipole_pos(:,3) ; dipole_pos(:,3) ; -dipole_pos(:,3) ; -dipole_pos(:,3) ];
% Concatenate decay rates for the whole 2pi angles
rad = [rad; rad ; rad ; rad];
tot = [tot; tot ; tot ; tot];
%% Final plot
%% plot 2d map of quantum yield enhancement
QY_x = rad(:,1)./(tot(:,1) + (1- Q_int)/Q_int ) ;
QY_y = rad(:,2)./(tot(:,2) + (1- Q_int)/Q_int ) ;
QY_z = rad(:,3)./(tot(:,3) + (1- Q_int)/Q_int ) ;
rad_average = (rad(:,1) + rad(:, 2) + rad(:, 3))./3 ;
tot_average = (tot(:,1) + tot(:, 2) + tot(: ,3))./3;
QY_average = rad_average./( tot_average + ( 1 - Q_int )/Q_int ); 
figure
subplot(2,2,1)
scatter3( pos(:,1), pos(:,2), QY_x, 8, QY_x, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\phi_\perp$'},'interpreter','latex')
subplot(2,2,2)
scatter3( pos(:,1), pos(:,2), QY_y, 8, QY_y, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\phi_{||}$'},'interpreter','latex')
subplot(2,2,3)
scatter3( pos(:,1), pos(:,2), QY_z, 8, QY_z, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\phi_{||}$'},'interpreter','latex')
subplot(2,2,4)
scatter3( pos(:,1), pos(:,2), QY_average, 8, QY_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\overline{\phi}$'},'interpreter','latex')
maximize(gcf)
saveas(gcf, 'QY enhancement.fig')
saveas(gcf, 'QY enhancement.png')

%% averaged decay rates
figure
subplot(1,2,1)
scatter3( pos(:,1), pos(:,2), rad_average, 8, rad_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\overline{\gamma^{rad}}/\gamma_0$'},'interpreter','latex')
subplot(1,2,2)
scatter3( pos(:,1), pos(:,2), tot_average, 8, tot_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\overline{\gamma^{tot}}/\gamma_0$'},'interpreter','latex')

saveas(gcf, 'decay_rates_averaged.fig')
saveas(gcf, 'decay_rates_averaged.png')
%% plot 2d map decay rates
figure
subplot(2,3,1)
scatter3( pos(:,1), pos(:,2), rad(: , 1 ),8, rad(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_\perp^{rad}/\gamma_0$'},'interpreter','latex')
subplot(2,3,2)
scatter3( pos(:,1), pos(:,2), rad(: , 2 ),8, rad(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')
subplot(2,3,3)
scatter3( pos(:,1), pos(:,2), rad(: , 3 ),8, rad(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')
subplot(2,3,4)
scatter3( pos(:,1), pos(:,2), tot(: , 1 ),8, tot(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_\perp^{tot}/\gamma_0$'},'interpreter','latex')
subplot(2,3,5)
scatter3( pos(:,1), pos(:,2), tot(: , 2 ),8, tot(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_{||}^{tot}/\gamma_0$'},'interpreter','latex')
subplot(2,3,6)
scatter3( pos(:,1), pos(:,2), tot(: , 3 ),8, tot(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([ 0 90 ])
title({'$\gamma_{||}^{tot}/\gamma_0$'},'interpreter','latex')
maximize(gcf)
saveas(gcf, 'decay_rates_components.fig')
saveas(gcf, 'decay_rates_components.png')
%% save data
save gamma_rad_BEM_2D rad
save gamma_tot_BEM_2D tot
save dipole_pos_2D dipole_pos
save pos_full pos
save QY_average_BEM_2D QY_average
save QY_x_2D QY_x
save QY_z_2D QY_z
save QY_y_2D QY_y

save particle p