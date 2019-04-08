%% BEM calculation of nanospheres

clear all
close all
clc

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  table of dielectric functions
epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };
%  diameter of sphere
diameter = 60;
Q_int = 0.01 ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 51 ),  ...
     pi * linspace( 0, 1, 51 ) .^ 2, diameter, 'triangles' );

% initialize particle
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
    
%%  dipole oscillator
enei = 550;

% 2D positions of dipole
[ x, z ] = SphereGrid2D ( diameter, 10, 50); 

pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 1.5 );

%  dipole excitation , x , y and z direction
dir = [ 1, 0, 0 ; 0 , 1 , 0 ; 0, 0, 1 ];


dip = dipole( pt, dir, op);

figure
bem_plot(p)
hold on
plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3) ,'r.')
view([0 0])
hold off
%%
saveas(gcf, 'dipole-sphere.fig')
saveas(gcf, 'dipole-sphere.png')
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
%%  final plot
%% plot 2d map of quantum yield enhancement
QY_x = rad(:,1)./(tot(:,1) + (1- Q_int)/Q_int ) ;
QY_y = rad(:,2)./(tot(:,2) + (1- Q_int)/Q_int ) ;
QY_z = rad(:,3)./(tot(:,3) + (1- Q_int)/Q_int ) ;
rad_average = (rad(:,1) + rad(:, 2) + rad(:, 3))./3 ;
tot_average = (tot(:,1) + tot(:, 2) + tot(: ,3))./3;
QY_average = (QY_x + QY_y + QY_z)./3 ;
%% 
figure
subplot(2,2,1)
scatter3( pos(:,1), pos(:,2), QY_x, 8, QY_x, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('x dipole orientation')
subplot(2,2,2)
scatter3( pos(:,1), pos(:,2), QY_y, 8, QY_y, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('y dipole orientation')
subplot(2,2,3)
scatter3( pos(:,1), pos(:,2), QY_z, 8, QY_z, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('z dipole orientation')
subplot(2,2,4)
scatter3( pos(:,1), pos(:,2), QY_average, 8, QY_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('average dipole orientation')
maximize(gcf)
saveas(gcf, 'QY enhancement_2D.fig')
saveas(gcf, 'QY enhancement_2D.png')
%% plot 2d map decay rates
figure
subplot(2,3,1)
scatter3( pos(:,1), pos(:,2), rad(: , 1 ),8, rad(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('x dipole radiative decay rate')
subplot(2,3,2)
scatter3( pos(:,1), pos(:,2), rad(: , 2 ),8, rad(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('y dipole radiative decay rate')
subplot(2,3,3)
scatter3( pos(:,1), pos(:,2), rad(: , 3 ),8, rad(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('z dipole radiative decay rate')
subplot(2,3,4)
scatter3( pos(:,1), pos(:,2), tot(: , 1 ),8, tot(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('x dipole total decay rate')
subplot(2,3,5)
scatter3( pos(:,1), pos(:,2), tot(: , 2 ),8, tot(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('y dipole total decay rate')
subplot(2,3,6)
scatter3( pos(:,1), pos(:,2), tot(: , 3 ),8, tot(:, 3),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title('z dipole total decay rate')
maximize(gcf)
saveas(gcf, 'decayrates_2D.fig')
saveas(gcf, 'decayrates_2D.png')
%% save data
save gamma_rad_BEM_2D rad
save gamma_tot_BEM_2D tot
save dipole_pos_2D dipole_pos
save dipole_pos_full_2D pos
save QY_average_BEM_2D QY_average
save QY_x_2D QY_x
save QY_z_2D QY_z
save QY_y_2D QY_y


%%  Post processing
%% QY maps
% Interpolate the 2D maps with finer grid and output the images
A_QY = scatteredInterpolant( pos( : ,1 ), pos(:, 2), QY_average );
A_QYx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), QY_x );
A_QYy = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), QY_y );
A_QYz = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), QY_z );
% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 400  ), linspace( -100, 100, 400 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 1.5 );
% Do interpolation and build interpolated results
QY_fined = A_QY( pt_fine.pos(:,1),pt_fine.pos(:,3) );
QYx_fined = A_QYx( pt_fine.pos(:,1),pt_fine.pos(:,3) );
QYy_fined = A_QYy( pt_fine.pos(:,1),pt_fine.pos(:,3) );
QYz_fined = A_QYz( pt_fine.pos(:,1),pt_fine.pos(:,3) );
%% 
figure
ax1 = subplot(2,2,1);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QY_fined(:),8,QY_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('QY Average dipole orientation')
ax2 = subplot(2,2,2);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYx_fined(:),8,QYx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('x dipole orientation')
ax3 = subplot(2,2,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYy_fined(:),8,QYy_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('y dipole orientation')
ax4 = subplot(2,2,4);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYz_fined(:),8,QYz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
title('z dipole orientation')

set(gcf,'position',[ 500, 200, 1000, 700 ])

saveas(gcf, 'QY_interp.fig')
saveas(gcf, 'QY_interp.png')
%% Decay rate maps
A_rad = scatteredInterpolant( pos( : ,1 ), pos(:, 2), rad_average );
A_radx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,1) );
A_rady = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,2) );
A_radz = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,3) );
A_tot = scatteredInterpolant( pos( : ,1 ), pos(:, 2), tot_average );
A_totx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,1) );
A_toty = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,2) );
A_totz = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,3) );
% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 200  ), linspace( -100, 100, 200 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 1.5 );
% Do interpolation and build interpolated results
rad_fined = A_rad( pt_fine.pos(:,1),pt_fine.pos(:,3) );
radx_fined = A_radx( pt_fine.pos(:,1),pt_fine.pos(:,3) );
rady_fined = A_rady( pt_fine.pos(:,1),pt_fine.pos(:,3) );
radz_fined = A_radz( pt_fine.pos(:,1),pt_fine.pos(:,3) );
tot_fined = A_tot( pt_fine.pos(:,1),pt_fine.pos(:,3) );
totx_fined = A_totx( pt_fine.pos(:,1),pt_fine.pos(:,3) );
toty_fined = A_toty( pt_fine.pos(:,1),pt_fine.pos(:,3) );
totz_fined = A_totz( pt_fine.pos(:,1),pt_fine.pos(:,3) );
%% 
figure
ax1 = subplot(2,4,1);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),rad_fined(:),8,rad_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Radiative decay rates average dipole orientation')
ax2 = subplot(2,4,2);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radx_fined(:),8,radx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Radiative decay rates  x dipole orientation')
ax3 = subplot(2,4,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),rady_fined(:),8,rady_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Radiative decay rates y dipole orientation')
ax4 = subplot(2,4,4);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radz_fined(:),8,radz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
title('Radiative decay rates z dipole orientation')

ax5 = subplot(2,4,5);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),tot_fined(:),8,tot_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('tot decay rates average dipole orientation')

ax6 = subplot(2,4,6);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totx_fined(:),8,totx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('tot decay rates x dipole orientation')

ax7 = subplot(2,4,7);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),toty_fined(:),8,toty_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('tot decay rates y dipole orientation')

ax8 = subplot(2,4,8);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totz_fined(:),8,totz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('tot decay rates z dipole orientation')
maximize
%% 
saveas(gcf, 'decayrates_interp.fig')
saveas(gcf, 'decayrates_interp.png')
save particle p
