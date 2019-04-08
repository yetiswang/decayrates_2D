%% Plot near field, quantum yield enhancement and total fluorescence enhancement map


Q_int = 0.01; 



%% average oriented dipole
figure
FE = ee_full.* QY_average/Q_int ;

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
epstab = { epsconst( 1.3^2 ), epstable( 'gold.dat' ) };
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
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
save FE_map FE


%% Plot total decay rate enhancement map
[row, col] = size(tot);
if col == 3
tot_avg = ( tot(:,1) + tot(:,2) + tot(:,3) ) ./3 ;
elseif col == 2
tot_avg = (tot(:,1) + 2.*tot(:, 2) )./3 ;
end 

figure
scatter3(pos(:,1),pos(:,2),tot_avg(:),5,tot_avg(:),'filled')
view([0 90])
xlabel('x (nm)')
ylabel('y (nm)')
title('Lifetime  ')
saveas(gcf,'lifetime.fig')
saveas(gcf,'lifetime.png')

%% Interpolate tot map for illustration
A_tot = scatteredInterpolant( pos( : ,1 ), pos(:, 2), tot_avg(:) );
mesh = [ 41, 41, 41]; % n1 for the circumference of the rod, n2 for the polar angles of the rod caps, n3 for the cylinder-shaped middle part of the rod
height =  70   ;
diameter =  30  ;
epstab = { epsconst( 1.3^2 ), epstable( 'gold.dat' ) };
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );
% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 600  ), linspace( -100, 100, 600 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
tot_fined = A_tot( pt_fine.pos(:,1),pt_fine.pos(:,3) );
%%
figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),tot_fined(:),8,tot_fined(:),'filled')
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
saveas(gcf,'lifetime_interp.fig')
saveas(gcf,'lifetime_interp.png')