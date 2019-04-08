%% Plotting results for BEM decay rates simulation
%  Input: 
%         rad: n*3 array, with the first dimension the perpendicular
%         dipoles
%         tot: n*3 array, with the first dimension the perpendicular
%         dipoles
%         pos: dipole positions
%         Q_int: intrinsic QY


function [ particle ] = DecayRatesPlot( p, op, pos, rad, tot, Q_int )

ndir = ['QY',num2str(Q_int)];
mkdir(ndir)
cd(ndir)

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
saveas(gcf, 'decay rates.fig')
saveas(gcf, 'decay rates.png')
%% save data
save gamma_rad_BEM_2D rad
save gamma_tot_BEM_2D tot
%save dipole_pos_2D dipole_pos
save pos_full pos
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
[ x_fine, z_fine ] = meshgrid( linspace( -150 , 150, 600  ), linspace( -150, 150, 600 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
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
title({'$\overline{\phi}$'},'interpreter','latex')
ax2 = subplot(2,2,2);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYx_fined(:),8,QYx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\phi_\perp$'},'interpreter','latex')
ax3 = subplot(2,2,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYy_fined(:),8,QYy_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\phi_{||}$'},'interpreter','latex')
ax4 = subplot(2,2,4);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYz_fined(:),8,QYz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
title({'$\phi_{||}$'},'interpreter','latex')

set(gcf,'position',[ 500, 200, 1000, 700 ])

saveas(gcf, 'QY_interp.fig')
saveas(gcf, 'QY_interp.png')
%% Decay rate maps
A_rad = scatteredInterpolant( pos( : ,1 ), pos(:, 2), rad_average, 'nearest' );
A_radx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,1), 'nearest' );
A_rady = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,2),'nearest' );
A_radz = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,3),'nearest' );
A_tot = scatteredInterpolant( pos( : ,1 ), pos(:, 2), tot_average,'nearest' );
A_totx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,1),'nearest' );
A_toty = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,2),'nearest' );
A_totz = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,3),'nearest' );
% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -150 , 150, 300  ), linspace( -150, 150, 300 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
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
title({'$\overline{\gamma}^{rad}/\gamma_0$'},'interpreter','latex')

ax2 = subplot(2,4,2);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radx_fined(:),8,radx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight

title({'$\gamma_\perp^{rad}/\gamma_0$'},'interpreter','latex')
ax3 = subplot(2,4,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),rady_fined(:),8,rady_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')
ax4 = subplot(2,4,4);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radz_fined(:),8,radz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')

ax5 = subplot(2,4,5);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),tot_fined(:),8,tot_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\overline{\gamma}^{tot}/\gamma_0$'},'interpreter','latex')

ax6 = subplot(2,4,6);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totx_fined(:),8,totx_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_\perp^{tot}/\gamma_0$'},'interpreter','latex')

ax7 = subplot(2,4,7);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),toty_fined(:),8,toty_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_{||}^{tot}/\gamma_0$'},'interpreter','latex')

ax8 = subplot(2,4,8);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totz_fined(:),8,totz_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_{||}^{tot}/\gamma_0$'},'interpreter','latex')
maximize
%% 
saveas(gcf, 'decayrates_interp.fig')
saveas(gcf, 'decayrates_interp.png')
save particle p