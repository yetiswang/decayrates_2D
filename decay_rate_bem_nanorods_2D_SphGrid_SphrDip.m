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
height =  42   ;
diameter =  30  ;

Q_int = [0.01, 0.05, 0.25, 0.75] ; % intrinsic quantum efficiency
%  nanorod with finer discretization at the top
%  To calculate decay rates close to surface, the mesh close to the
%  positions of dipoles need to refined.
p = trirod ( diameter, height, mesh, 'triangles' );
%  initialize sphere
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
%  rotate the particle for 90 degrees
p = rot (p, 90, [0, -1, 0]);
%%  dipole oscillator
% 2D positions of dipole
% Consuming large amount of memory, use mesh size in the units of nanometer
% inbed points on the side of the nanorods
[ x, z ] = RodGrid2D ( diameter, height, 60, 50 );
pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
% calculate the spherical dipole orientation for each dipole location
dir = GetDipDirGNR( height, diameter, pt );
dip = dipoleret( pt, dir, 'full');

%% Build new direcotry
parent = pwd;
directory = ['H',num2str(height),'D',num2str(diameter),'_fluorescence_enhancement'];
mkdir(directory)
cd(directory)
[ sca, fit, Lorentz ] = spect_GNR_BEM( epstab, height, diameter, linspace(500,800,20 ));
enei = round(1248./Lorentz(3)) + 25 ; % 25 nm red to the resonance wavelength
vline(enei, 'r-','\lambda_{dip}');
saveas (gcf, ['spectrum.fig'], 'fig')
saveas (gcf, ['spectrum.png'], 'png')
sub_dir = ['Dip',num2str(enei),'nm'];
mkdir(sub_dir)
cd(sub_dir)
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

for i = 1 : length (Q_int)
    cd(parent)
    cd(directory)
    cd(sub_dir)
    DecayRatesPlot( p, op, pos, rad, tot, Q_int(i) );
end

%% calculation of near field excitation enhancement 
cd(parent)
cd(directory)
[ee_full] = GNR_NF_2D( epstab, height, diameter, mesh, pt, op ); 

%% Hallelujah
finishedScript