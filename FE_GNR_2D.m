%% BEM dipole decay rate and quantum yield calculation in the vinicity of nanorods

function FE_GNR_2D( metal, height, diameter, mesh, scale, enei_field, enei_dipole )

%%  initialization
%  options for BEM simulation
op = bemoptions( 'sim', 'ret', 'interp', 'curv' );

switch metal
    case 'gold'
        %  table of dielectric functions
        epstab = { epsconst( 1.33^2 ), epstable( 'gold.dat' ) };
        if exist('enei_field','var') == 0 || isempty(enei_field)
            directory = ['H',num2str(height),'D',num2str(diameter),'_Au_fluorescence_enhancement'];
        else
            directory = ['H',num2str(height),'D',num2str(diameter),'_Au_fluorescence_enhancement','_Exc',num2str(enei_field),'_Dip', num2str(enei_dipole)];
        end
    case 'silver'
        epstab = { epsconst( 1.33^2 ), epstable( 'silver.dat' ) };
        if exist('enei_field','var') == 0 || isempty(enei_field)
            directory = ['H',num2str(height),'D',num2str(diameter),'_Ag_fluorescence_enhancement'];
        else
            directory = ['H',num2str(height),'D',num2str(diameter),'_Ag_fluorescence_enhancement','_Exc',num2str(enei_field),'_Dip', num2str(enei_dipole)];
        end
end

Q_int = [0.05, 0.25, 0.4, 0.75] ; % intrinsic quantum efficiency
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
[ x, z ] = RodGrid2D ( diameter, height, 60, 20, scale );
pt = compoint( p, [ x( : ), 0 * x( : ) , z( : ) ], op, 'medium', 1, 'mindist' , 0.8 );
% calculate the spherical dipole orientation for each dipole location
dir = GetDipDirGNR( height, diameter, pt );
dip = dipoleret( pt, dir, 'full');

%% Build new direcotry
parent = pwd;
mkdir(directory)
cd(directory)
[ sca, fit, Lorentz ] = spect_GNR_BEM( epstab, height, diameter, linspace(300,850,40 ));
if exist('enei_dipole','var') == 0 || isempty(enei_dipole)
    enei_dipole = round(1248./Lorentz(3)) + 25 ; % 25 nm red to the resonance wavelength
else
end
vline(enei_dipole, 'r-','\lambda_{dip}');
saveas (gcf, ['spectrum.fig'], 'fig')
saveas (gcf, ['spectrum.png'], 'png')
sub_dir = ['Dip',num2str(enei_dipole),'nm'];
mkdir(sub_dir)
cd(sub_dir)
%%
figure
plot(p)
hold on
coneplot(pt.pos,dir)
view([0 0])
axis on
grid on
hold off
%%
saveas(gcf,'dipole-orientation.fig')
saveas(gcf,'dipole-orientation.png')
%%
figure
plot(p)
hold on
plot3( pt.pos(:,1),pt.pos(:,2),pt.pos(:,3),'r*')
hold off
view([0 0])

%%
saveas(gcf,'dipole-sphere.fig')
saveas(gcf,'dipole-sphere.png')
%%  BEM simulation
%  set up BEM solver
bem = bemsolver( p, op );
%  surface charge
sig = bem \ dip( p, enei_dipole );
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
close all
%% calculation of near field excitation enhancement
cd(parent)
cd(directory)

if exist('enei_field','var') == 0 || isempty(enei_field)
    [ee_full] = GNR_NF_2D( epstab, height, diameter, mesh, pt, op ); %% choose excitation wavelength
else
    [ee_full] = GNR_NF_2D( epstab, height, diameter, mesh, pt, op, x, z, enei_field ); %% choose excitation wavelength
end

close all
%% Hallelujah
finishedScript
end