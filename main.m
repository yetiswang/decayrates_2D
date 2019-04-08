%% main script for running GNR fluorescence

function main

clear all
close all
height = [ 60 ];
diameter = [ 30 ];

parent = pwd;
for i = 1 : length(height)
    if height(i) < 70 && ~isequal(height(i),diameter(i))
        mesh = [ 41 41 41 ];
        scale = 1;
        particle = 'A nanoroooood!';
        disp(particle);
    elseif height(i) >= 70 && ~isequal(height(i),diameter(i))
        mesh = [ 45 45 60 ];
        scale = 0.3 ;
        particle = 'A nanoroooood!';
        disp(particle);
    elseif isequal(height(i),diameter(i))
        particle = 'A nanospheeeeeere!';
        disp(particle);
    end
    
    switch particle
        case 'A nanoroooood!'
            FE_GNR_2D( 'gold', height(i), diameter(i), mesh, scale, 637, 670 )
        case  'A nanospheeeeeere!'
            FE_GNS_2D( diameter(i) )
    end
    cd(parent)
    close all
end

end
