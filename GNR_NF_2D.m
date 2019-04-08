%% Calculate the near field
%  This script work with decay rates calculation of gold nanorods.
%  Given the dielectric function, shape paramters, mesh size and dipole positions
%  exactly in the same manner of decay rates calculation, the near field enhancement
%  can be calculated and saved.

function [ee_full] = GNR_NF_2D( epstab, height, diameter, mesh, pt, op, x, z, enei_field )

ndir = ['Exc_H',num2str(height),'_D',num2str(diameter) ];
mkdir(ndir)
cd(ndir)

%% LOOP FOR PARTICLE GEOMETRIES
%  initialize nanosphere
p = trirod ( diameter, height, mesh, 'triangles' );

%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);

%  plot particle boundary
figure
plot( p, 'EdgeColor', 'b');
axis tight;
grid on;

saveas (gcf, ['geometry.fig'], 'fig')
saveas (gcf, ['geometry.png'], 'png')
close gcf

%  set up COMPARTICLE object
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );

%  TM mode, excitation from above
dir = [ 0, 0, -1];
pol = [ 1, 0, 0 ];

%% NEAR FIELD ENHANCEMENT
[ sca, fit, Lorentz ] = spect_GNR_BEM( epstab, height, diameter, linspace(400,1000,50 ));
if nargin == 6
    enei_field = round(1248./Lorentz(3)) ; %  resonance wavelength
else
end
vline(enei_field, 'r-','\lambda_{exc}');
saveas (gcf, ['spectrum.fig'], 'fig')
saveas (gcf, ['spectrum.png'], 'png')
close gcf
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
%%  computation of electric field
%  object for electric field
%    MINDIST controls the minimal distance of the field points to the
%    particle boundary, MESHFIELD must receive the OP structure which also
%    stores the table of precomputed reflected Green functions
emesh = meshfield( p, x( : ), 0 * x( : ) , z( : ), op, 'mindist', 0.2, 'nmax', 2000 , 'waibar', 1 );
%  induced and incoming electric field
e = emesh( sig ) + emesh( exc.field( emesh.pt, enei_field ) );
%  norm of electric field
enorm = vecnorm(e);

ee = enorm.^2;
%% plot near field

%  get electric field
ex = sqrt( dot (e(: ,1) , e(: ,1), 	2 ) );
ey = sqrt( dot (e(: ,2) , e(: ,2), 	2 ) );
ez = sqrt( dot (e(: ,3) , e(: ,3), 	2 ) );

% plot electric field vector
figure
coneplot( pt.pos, e )
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
ylabel('z (nm)')
view([ 0 90 ])

title('Near field intensity')
saveas(gcf, 'nearfield_2D_full.fig')
saveas(gcf, 'nearfield_2D_full.png')

figure
subplot( 1, 3, 1 )
scatter3(pos(:,1),pos(:,2),ex_full(:),5,ex_full(:),'filled')
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('z (nm)')
legend('E_x^2/E_0^2')
view([ 0 90 ])

subplot( 1 ,3 ,2 )
scatter3(pos(:,1),pos(:,2),ey_full(:),5,ey_full(:),'filled')
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('z (nm)')
legend('E_y^2/E_0^2')
view([ 0 90 ])

subplot( 1 ,3 ,3 )
scatter3(pos(:,1),pos(:,2),ez_full(:),5,ez_full(:),'filled')
colorbar;
colormap jet( 1000 );
xlabel('x (nm)')
ylabel('z (nm)')
legend('E_z^2/E_0^2')
view([ 0 90 ])

title('Near field intensity')
set(gcf,'position',[99 375 1779 387])
saveas(gcf, 'nearfield_2Dxyz_full.fig')
saveas(gcf, 'nearfield_2Dxyz_full.png')

%% Interpolation for finer meshes for average full electric field enhancement, and electric field of each component
A = scatteredInterpolant( pos( : ,1 ), pos(:, 2), ee_full(:) );
% Build fined mesh grid for good looking pictures

[ x_fine, z_fine ] = meshgrid( linspace( -150 , 150, 600  ), linspace( -150, 150, 600 ) );
pt_fine = compoint( p, [ x_fine( : ), 0 * x_fine( : ) , z_fine( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
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
legend('E_x^2/E_0^2')
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
legend('E_y^2/E_0^2')
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
legend('E_z^2/E_0^2')
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


