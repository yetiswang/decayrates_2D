%% BEM calculation of nanospheres

function FE_GNS_2D( metal, diameter )
%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

switch metal
    case 'gold'
        %  table of dielectric functions
        epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };
        directory = ['D',num2str(diameter),'_Au_fluorescence_enhancement'];
    case 'silver'
        epstab = { epsconst( 1.33^2 ), epstable( 'silver.dat' ) };
        directory = ['D',num2str(diameter),'_Ag_fluorescence_enhancement'];
end

%  table of dielectric functions
epstab = { epsconst( 1.33^2 ), epstable( 'silver.dat' ) };

Q_int = [0.01, 0.05, 0.25, 0.75] ; % intrinsic quantum efficiency
%  nanosphere with finer discretization at the top
%  To calculate decay rates close to spheres, the mesh close to the
%  positions of dipoles need to refined.
p = trispheresegment( 2 * pi * linspace( 0, 1, 51 ),...
    pi * linspace( 0, 1, 51 ) .^ 2, diameter, 'triangles' );

% initialize particle
p = comparticle( epstab, { p }, [ 2, 1 ], 1, op );
%% make a new directory
parent = pwd;
mkdir(directory)
cd(directory)
%% Calculate and save the scattering spectrum for reference
[ sca, fit, Lorentz ] = spect_GNS_BEM( epstab, diameter, linspace(300,800,30 ));
%%  dipole oscillator
enei = round(1248./Lorentz(3)) + 25 ; % 25 nm red to the resonance wavelength
vline(enei, 'r-','\lambda_{dip}');
saveas (gcf, ['spectrum.fig'], 'fig')
saveas (gcf, ['spectrum.png'], 'png')
sub_dir = ['Dip',num2str(enei),'nm'];
mkdir(sub_dir)
cd(sub_dir)

% 2D positions of dipole
[ x, z ] = SphereGrid2D ( diameter, 50, 50);

pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 1.5 );

%  get radial and tangential orientation of the dipole for each dipole
%  position

for i = 1 : numel(pt.pos(:,1))
    vecnorm = sqrt ( pt.pos(i,1)^2 + pt.pos(i,2)^2 + pt.pos(i,3)^2 ) ;
    dir(i,1,1) = - pt.pos(i,1)/vecnorm;
    dir(i,1,2) = 0;
    dir(i,2,1) = - pt.pos(i,2)/vecnorm;
    dir(i,2,2) = 1;
    dir(i,3,1) = - pt.pos(i,3)/vecnorm;
    dir(i,3,2) = 0;
end

dip = dipoleret( pt, dir, 'full');

figure
bem_plot(p)
hold on
plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3) ,'r.')
view([0 0])
hold off
%%
saveas(gcf, 'dipole-sphere.fig')
saveas(gcf, 'dipole-sphere.png')
%%
figure
coneplot(pt.pos,dir)
axis on
grid on
%%
saveas(gcf,'dipole-orientation.fig')
saveas(gcf,'dipole-orientation.png')

%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei );
%  total and radiative decay rate
[ tot, rad ] = dip.decayrate( sig );

%% Expand to full angle sphere slice
% Concatenate position of dipoles for the whole 2pi angles
dipole_pos = pt.pos ;
pos(:,1) = [ dipole_pos(:,1) ; ( - dipole_pos(:,1) ) ; ( - dipole_pos(:,1) ) ; dipole_pos(:,1)];
pos(:,2) = [dipole_pos(:,3) ; dipole_pos(:,3) ; -dipole_pos(:,3) ; -dipole_pos(:,3) ];
% Concatenate decay rates for the whole 2pi angles
rad = [rad; rad ; rad ; rad];
tot = [tot; tot ; tot ; tot];

sub_dir = ['Dip',num2str(enei),'nm'];
mkdir(sub_dir)
cd(sub_dir)
%% Final plot

for i = 1 : length (Q_int)
    cd(parent)
    cd(directory)
    cd(sub_dir)
    SphereDecayRatesPlot( p, op, pos, rad, tot, Q_int(i) );
end
%% calculation of near field excitation enhancement 
cd(parent)
cd(directory)
[ee_full] = GNS_NF_2D( epstab, diameter, pt, op );

%% Hallelujah
finishedScript
close all
end 
