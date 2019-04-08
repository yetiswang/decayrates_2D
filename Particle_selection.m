%% Selection of gold nanorods based on width, length of SPR wavelength

clear all
close all
%%
diameter = [ 10 10 10 10 60 60 60 60  ];
height = [  20 25 28 30 125 127 129 130];
epstab = { epsconst( 1.33^2 ), epstable( 'silver.dat' ) };

parent = pwd;
for i = 1 : length(height)
    if height(i) == diameter(i)
        ndir = ['D',num2str(diameter(i))];
        mkdir(ndir)
        cd(ndir)
        [ sca, fit, Lorentz ] = spect_GNS_BEM( epstab, diameter(i), linspace(300,700,30 ));
        cd(parent)
    else
        
        ndir = ['H',num2str(height(i)),'D',num2str(diameter(i))];
        mkdir(ndir)
        cd(ndir)
        [ sca, fit, Lorentz ] = spect_GNR_BEM( epstab, height(i), diameter(i), linspace(400,1000,30 ));
        cd(parent)
    end
end