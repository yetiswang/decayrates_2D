%% plot decay rates of sphere-dipoles

function SphereDecayRatesPlot( p, op, pos, rad, tot, Q_int )

ndir = ['QY',num2str(Q_int)];
mkdir(ndir)
cd(ndir)
%%  final plot
%% plot 2d map of quantum yield enhancement
QY_perp = rad(:,1)./(tot(:,1) + (1- Q_int)/Q_int ) ;
QY_para = rad(:,2)./(tot(:,2) + (1- Q_int)/Q_int ) ;
rad_average = (rad(:,1) + 2.*rad(:, 2) )./3 ;
tot_average = (tot(:,1) + 2.*tot(:, 2) )./3 ;
QY_average = rad_average./( tot_average + (1- Q_int)/Q_int  ) ;
%% 
figure
subplot(1,3,1)
scatter3( pos(:,1), pos(:,2), QY_perp, 8, QY_perp, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\phi_\perp$'},'interpreter','latex')
subplot(1,3,2)
scatter3( pos(:,1), pos(:,2), QY_para, 8, QY_para, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\phi_{||}$'},'interpreter','latex')
subplot(1,3,3)
scatter3( pos(:,1), pos(:,2), QY_average, 8, QY_average, 'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\overline{\phi}$'},'interpreter','latex')
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
title({'$\gamma_\perp^{rad}/\gamma_0$'},'interpreter','latex')
subplot(2,3,2)
scatter3( pos(:,1), pos(:,2), rad(: , 2 ),8, rad(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')
subplot(2,3,3)
scatter3( pos(:,1), pos(:,2), tot(: , 1 ),8, tot(:, 1),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\gamma_\perp^{tot}/\gamma_0$'},'interpreter','latex')
subplot(2,3,4)
scatter3( pos(:,1), pos(:,2), tot(: , 2 ),8, tot(:, 2),'filled' )
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
view([0 90])
title({'$\gamma_{||}^{tot}/\gamma_0$'},'interpreter','latex')
maximize(gcf)
saveas(gcf, 'decayrates_2D.fig')
saveas(gcf, 'decayrates_2D.png')
%% save data
save gamma_rad_BEM_2D rad
save gamma_tot_BEM_2D tot
%save dipole_pos_2D dipole_pos
save dipole_pos_full_2D pos
save QY_average_BEM_2D QY_average
save QY_perp_2D QY_perp
save QY_para_2D QY_para

%%  Post processing
%% QY maps
% Interpolate the 2D maps with finer grid and output the images
A_QY = scatteredInterpolant( pos( : ,1 ), pos(:, 2), QY_average );
A_QYperp = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), QY_perp );
A_QYpara = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), QY_para );

% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 400  ), linspace( -100, 100, 400 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 1.5 );
% Do interpolation and build interpolated results
QY_fined = A_QY( pt_fine.pos(:,1),pt_fine.pos(:,3) );
QYperp_fined = A_QYperp( pt_fine.pos(:,1),pt_fine.pos(:,3) );
QYpara_fined = A_QYpara( pt_fine.pos(:,1),pt_fine.pos(:,3) );

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
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYperp_fined(:),8,QYperp_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\phi_\perp$'},'interpreter','latex')
ax3 = subplot(2,2,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),QYpara_fined(:),8,QYpara_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\phi_{||}$'},'interpreter','latex')

set(gcf,'position',[ 500, 200, 1000, 700 ])

saveas(gcf, 'QY_interp.fig')
saveas(gcf, 'QY_interp.png')
%% Decay rate maps
A_rad = scatteredInterpolant( pos( : ,1 ), pos(:, 2), rad_average );
A_radperp = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,1) );
A_radpara = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), rad(:,2) );

A_tot = scatteredInterpolant( pos( : ,1 ), pos(:, 2), tot_average );
A_totx = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,1) );
A_toty = scatteredInterpolant( pos( :, 1 ), pos(:, 2 ), tot(:,2) );

% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 200  ), linspace( -100, 100, 200 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 1.5 );
% Do interpolation and build interpolated results
rad_fined = A_rad( pt_fine.pos(:,1),pt_fine.pos(:,3) );
radperp_fined = A_radperp( pt_fine.pos(:,1),pt_fine.pos(:,3) );
radpara_fined = A_radpara( pt_fine.pos(:,1),pt_fine.pos(:,3) );

tot_fined = A_tot( pt_fine.pos(:,1),pt_fine.pos(:,3) );
totperp_fined = A_totx( pt_fine.pos(:,1),pt_fine.pos(:,3) );
totpara_fined = A_toty( pt_fine.pos(:,1),pt_fine.pos(:,3) );

%% 
figure
ax1 = subplot(2,3,1);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),rad_fined(:),8,rad_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\overline{\gamma}^{rad}/\gamma_0$'},'interpreter','latex')
ax2 = subplot(2,3,2);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radperp_fined(:),8,radperp_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_\perp^{rad}/\gamma_0$'},'interpreter','latex')
ax3 = subplot(2,3,3);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),radpara_fined(:),8,radpara_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_{||}^{rad}/\gamma_0$'},'interpreter','latex')
ax5 = subplot(2,3,4);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),tot_fined(:),8,tot_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\overline{\gamma}^{tot}/\gamma_0$'},'interpreter','latex')

ax6 = subplot(2,3,5);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totperp_fined(:),8,totperp_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title({'$\gamma_\perp^{tot}/\gamma_0$'},'interpreter','latex')

ax7 = subplot(2,3,6);
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),totpara_fined(:),8,totpara_fined(:),'filled')
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