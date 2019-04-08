%% Plot near field, quantum yield enhancement and total fluorescence enhancement map


Q_int = 0.01; 

%% X oriented dipole
figure
subplot(1,3,1)
scatter3( pos(:,1), pos(:,2), ex_full(:), 5, ex_full(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity (x-oriented dipole)')
subplot(1,3,2)
scatter3( pos(:,1), pos(:,2), QY_x(:), 5, QY_x(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('QY enhancement (x-oriented dipole)')
% Combine NF and QY enhancement
FE_x = ex_full .* ( QY_x./Q_int ) ; 
subplot(1,3,3)
scatter3( pos(:,1), pos(:,2), FE_x(:), 5, FE_x(:), 'filled')
colorbar
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
view( [0 90] ) 
title('Fluorescence enhancement (x-oriented dipole)')
set(gcf,'position', [ 62 600 1808 379 ])

saveas(gcf, 'fluorescence_enhancement_x.fig')
saveas(gcf, 'fluorescence_enhancement_x.png')
%% Y oriented dipole
figure
subplot(1,3,1)
scatter3( pos(:,1), pos(:,2), ey_full(:), 5, ey_full(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity (y-oriented dipole)')
subplot(1,3,2)
scatter3( pos(:,1), pos(:,2), QY_y(:), 5, QY_y(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('QY enhancement (y-oriented dipole)')
% Combine NF and QY enhancement
FE_y = ey_full .* ( QY_y./Q_int ) ; 
subplot(1,3,3)
scatter3( pos(:,1), pos(:,2), FE_y(:), 5, FE_y(:), 'filled')
colorbar
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
view( [0 90] ) 
title('Fluorescence enhancement (y-oriented dipole)')
set(gcf,'position', [ 62 600 1808 379 ])

saveas(gcf, 'fluorescence_enhancement_y.fig')
saveas(gcf, 'fluorescence_enhancement_y.png')
%% z oriented dipole
figure
subplot(1,3,1)
scatter3( pos(:,1), pos(:,2), ez_full(:), 5, ez_full(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity (z-oriented dipole)')
subplot(1,3,2)
scatter3( pos(:,1), pos(:,2), QY_z(:), 5, QY_z(:), 'filled')
colorbar
colormap jet( 1000 );
view( [0 90] ) 
xlabel('x (nm)')
ylabel('y (nm)')
title('QY enhancement (z-oriented dipole)')
% Combine NF and QY enhancement
FE_z = ez_full .* ( QY_z./Q_int ) ; 
subplot(1,3,3)
scatter3( pos(:,1), pos(:,2), FE_z(:), 5, FE_z(:), 'filled')
colorbar
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
view( [0 90] ) 
title('Fluorescence enhancement (z-oriented dipole)')
set(gcf,'position', [ 62 600 1808 379 ])

saveas(gcf, 'fluorescence_enhancement_z.fig')
saveas(gcf, 'fluorescence_enhancement_z.png')
%% average oriented dipole
figure
FE = ( FE_x + FE_y + FE_z )./3 ;
scatter3( pos(:,1), pos(:,2), FE(:), 5, FE(:), 'filled')
colorbar
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
view( [0 90] ) 
title('Fluorescence enhancement (average-oriented dipole)')

saveas(gcf, 'fluorescence_enhancement.fig')
saveas(gcf, 'fluorescence_enhancement.png')

%% Interpolate FE map for illustration
A_EF = scatteredInterpolant( pos( : ,1 ), pos(:, 2), FE(:) );
mesh = [ 41, 41, 41]; % n1 for the circumference of the rod, n2 for the polar angles of the rod caps, n3 for the cylinder-shaped middle part of the rod
height =  70   ;
diameter =  30  ;
epstab = { epsconst( 1.3^2 ), epstable( 'gold.dat' ) };op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 600  ), linspace( -100, 100, 600 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
EF_fined = A_EF( pt_fine.pos(:,1),pt_fine.pos(:,3) );
%%
figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),EF_fined(:),8,EF_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
hold on 
plot(p)
hold off
title('Fluorescence enhancement average dipole orientation')
saveas(gcf,'FE_map_interp.fig')
saveas(gcf,'FE_map_interp.png')
%% 
save dipole_pos_full pos
save FE_map_x FE_x
save FE_map_y FE_y
save FE_map_z FE_z
save FE_map FE


%% Plot total decay rate enhancement map
tot_avg = ( tot(:,1) + tot(:,2) + tot(:,3) ) ./3 ;
scatter3(pos(:,1),pos(:,2),tot_avg(:),5,tot_avg(:),'filled')
view([0 90])
xlabel('x (nm)')
ylabel('y (nm)')
title('Lifetime  ')