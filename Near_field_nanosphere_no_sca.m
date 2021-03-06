%% Calculate the near field

%% clear workspace
clearvars
close all
clc
dbstop if error

%% initialization
eps_m = 1.33^2 ;
%  table of dielectric functions
epstab = { epsconst( eps_m ), epstable( 'gold.dat' ) };

%  options for BEM simulations
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

%  diameter of sphere
diameter = 30;
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 51 ),  ...
    pi * linspace( 0, 1, 51 ) .^ 2, diameter, 'triangles' );

% initialize particle
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  plot particle boundary
figure
plot( p, 'EdgeColor', 'b');
axis tight;
grid on;

%  TM mode, excitation from above
dir = [ 0, 0, -1];
pol = [ 1, 0, 0 ];

%% NEAR FIELD ENHANCEMENT
%% Calculate and save the scattering spectrum for reference
[ sca, fit, Lorentz ] = spect_GNS_BEM( epstab, diameter, linspace(450,650,20 ));
enei_field = round(1248./Lorentz(3)) ; % 25 nm red to the resonance wavelength
vline(enei_field, 'r-','\lambda_{exc}');
saveas (gcf, ['spectrum.fig'], 'fig')
saveas (gcf, ['spectrum.png'], 'png')

%%   BEM solver
%  initialize BEM solver
bem = bemsolver( p, op );
%  initialize plane wave excitation
exc = planewave( pol, dir, op );
%  solve BEM equation
sig = bem \ exc( p, enei_field );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );
%  close waitbar
multiWaitbar( 'CloseAll' );

% 2D positions of dipole
[ x, z ] = SphereGrid2D ( diameter, 50, 50);
pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 1.5 );

%  set up Green function object between PT and P
g = greenfunction( pt, p, op );
%  compute electric field
f = field( g, sig );

%% plot near field
%  plot electric field

ee = sqrt( dot (f.e , f.e, 2 ) ) ;
ee = ee.^2 ;

%  get electric field
ex = sqrt( dot (f.e(: ,1) , f.e(: ,1), 	2 ) );
ey = sqrt( dot (f.e(: ,2) , f.e(: ,2), 	2 ) );
ez = sqrt( dot (f.e(: ,3) , f.e(: ,3), 	2 ) );

% plot electric field vector
figure
coneplot( pt.pos, f.e )
axis on
grid on
hold on;
plot(p)
saveas( gcf, 'Electric field vector.fig' )
saveas( gcf, 'Electric field vector.png' )

% plot enhanced field in every component
figure
subplot(1,3,1)
scatter3(pt.pos(:,1),pt.pos(:,3),(ex(:)).^2, 5, (ex(:)).^2, 'filled'   );

colorbar;
colormap jet( 1000 );
view([0 90])
legend('E_x^2/E_0^2')
xlabel('x (nm)')
ylabel('z (nm)')


subplot(1,3,2)
scatter3(pt.pos(:,1),pt.pos(:,3),(ey(:)).^2, 5, (ey(:)).^2, 'filled'   );
colorbar;
colormap jet( 1000 );
view([0 90])
legend('E_y^2/E_0^2')
xlabel('x (nm)')
ylabel('z (nm)')

subplot(1,3,3)
scatter3(pt.pos(:,1),pt.pos(:,3),(ez(:)).^2, 5, (ez(:)).^2, 'filled'   );
colorbar;
colormap jet( 1000 );
view([0 90])
legend('E_z^2/E_0^2')
xlabel('x (nm)')
ylabel('z (nm)')
maximize

saveas( gcf, 'nearfield_xyz.fig' )
saveas( gcf, 'nearfield_xyz.png' )

figure
scatter3(pt.pos(:,1),pt.pos(:,3),ee(:),5,ee(:),'filled')
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
view( [ 0 90 ] )
title('Near field intensity')
saveas(gcf, 'nearfield_2D_quarter.fig')
saveas(gcf, 'nearfield_2D_quarter.png')


%% Expand the quarter sphere to full sphere
% Concatenate position of dipoles for the whole 2pi angles
% The postions for near field calculation is also the positions for
% dipoles in decay rates calculation

dipole_pos = pt.pos ;
pos(:,1) = [ dipole_pos(:,1) ; ( - dipole_pos(:,1) ) ; ( - dipole_pos(:,1) ) ; dipole_pos(:,1)];
pos(:,2) = [dipole_pos(:,3) ; dipole_pos(:,3) ; -dipole_pos(:,3) ; -dipole_pos(:,3) ];
ee_full = [ ee; ee; ee; ee ];
ex_full = [ ex; ex; ex; ex ];
ex_full = ex_full.^2;
ey_full = [ ey; ey; ey; ey ];
ey_full = ey_full.^2;
ez_full = [ ez; ez; ez; ez ];
ez_full = ez_full.^2;

%% plot full near field

figure
scatter3(pos(:,1),pos(:,2),ee_full(:),5,ee_full(:),'filled')
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('y (nm)')
title('Near field intensity')
saveas(gcf, 'nearfield_2D_full.fig')
saveas(gcf, 'nearfield_2D_full.png')

%% Interpolation for finer meshes
A = scatteredInterpolant( pos( : ,1 ), pos(:, 2), ee_full(:) );
% Build fined mesh grid for good looking pictures
[ x_fine, z_fine ] = meshgrid( linspace( -100 , 100, 400  ), linspace( -100, 100, 400 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 1.5 );
ee_fined = A( pt_fine.pos(:,1),pt_fine.pos(:,3) );

figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),ee_fined(:),8,ee_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Electric field')
saveas(gcf, 'nearfield_2D_full_interp.fig')
saveas(gcf, 'nearfield_2D_full_interp.png')

%% X component
A_x = scatteredInterpolant( pos( : ,1 ), pos(:, 2), ex_full(:) );
ex_fined = A_x( pt_fine.pos(:,1),pt_fine.pos(:,3) );
figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),ex_fined(:),8,ex_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Electric field')
saveas(gcf, 'nearfield_2Dx_full_interp.fig')
saveas(gcf, 'nearfield_2Dx_full_interp.png')

%% Y component
A_y = scatteredInterpolant( pos( : ,1 ), pos(:, 2), ey_full(:) );
ey_fined = A_y( pt_fine.pos(:,1),pt_fine.pos(:,3) );
figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),ey_fined(:),8,ey_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Electric field')
saveas(gcf, 'nearfield_2Dy_full_interp.fig')
saveas(gcf, 'nearfield_2Dy_full_interp.png')

%% Z component
A_z = scatteredInterpolant( pos( : ,1 ), pos(:, 2), ez_full(:) );
ez_fined = A_z( pt_fine.pos(:,1),pt_fine.pos(:,3) );
figure
scatter3(pt_fine.pos(:,1),pt_fine.pos(:,3),ez_fined(:),8,ez_fined(:),'filled')
view([0 90])
colormap jet
colorbar
xlabel('x (nm)')
ylabel('z (nm)')
axis tight
title('Electric field')
saveas(gcf, 'nearfield_2Dz_full_interp.fig')
saveas(gcf, 'nearfield_2Dz_full_interp.png')


%% Save coordinates and near field intensities
% Coordinates here should have the same coordinates as in the calculation
% for decay rates.
dipole_pos = pt.pos ;
save nearfieldquarter ee
save nearfieldx ex_full
save nearfieldy ey_full
save nearfeildz ez_full
save dipole_pos dipole_pos
save nearfieldfull ee_full
save dipole_pos_full_2D pos
