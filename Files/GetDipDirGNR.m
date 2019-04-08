%% Given a nanorod geometry, diameter and height, and dipole locations, calculate the 
%  the radial and tangential orientation of the dipole for each dipole
%  location

function dir = GetDipDirGNR( height,diameter,pt )
for i = 1 : numel(pt.pos(:,1))
    if pt.pos(i,1) <= (height/2 - diameter/2)
        %if the coordinate is outside of the cylinder part of the rod, and
        %on the positive side
        dir(i,1,1) = 0;
        dir(i,1,2) = 0;
        dir(i,1,3) = - 1;
        dir(i,2,1) = 0;
        dir(i,2,2) = 1;
        dir(i,2,3) = 0;
        dir(i,3,1) = - 1;
        dir(i,3,2) = 0;
        dir(i,3,3) = 0;
    elseif pt.pos(i,1) > (height/2 - diameter/2)
        vecnorm = sqrt ( (pt.pos(i,1) - (height/2 - diameter/2))^2 + pt.pos(i,2)^2 + pt.pos(i,3)^2 ) ;
        dir(i,1,1) = - ( pt.pos(i,1) - (height/2 - diameter/2))/vecnorm;
        dir(i,2,1) = - pt.pos(i,2)/vecnorm;
        dir(i,3,1) = - pt.pos(i,3)/vecnorm;
        perp = null(dir(i,:,1)).'; % get the perpendicular vectors
        dir(i,1,2) = perp(1,1);
        dir(i,2,2) = perp(1,2);
        dir(i,3,2) = perp(1,3);
        dir(i,1,3) = perp(2,1);
        dir(i,2,3) = perp(2,2);
        dir(i,3,3) = perp(2,3);
    else
    end
end
end 