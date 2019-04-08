%% Plots for finely interpolated maps as used in the ACS photonics paper
% Before running, import nearfieldfull, QY_average, and define shape and
% size of the particle and the QY0

height = 65 ;
diameter = 30 ;
QY_0 = 0.01; 
ExcludedParticle = [ diameter/2, height - diameter, 1 ];
%% Confocal + enhanced near field plot

A_exc = scatteredInterpolant( pos(:, 1), pos(:, 2), ee_full(:) );

[xx, yy] = meshgrid(linspace( -600e-9,600e-9,2000 ),linspace( -1000e-9,1000e-9, 2000) );
I_exc = zeros(size(xx));
%I_gauss = zeros(size(xx));

for i = 1 : numel(xx)
    I_exc(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_exc,ExcludedParticle );
end
% for i = 1 : numel(xx)
%      I_gauss(i) = D3GaussFunction_analytical( confocal.params,[xx(i),0,yy(i)]);
% end

figure
contours = [ 1 10 100 1000 ];
contour(I_exc);
% hold on
% contour(I_gauss);

figure
%I_tot = I_exc + I_gauss ;
I_tot = I_exc;
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],log10(I_tot))
colormap jet
hold on
y = get(colorbar,'YTick');
colorbar('YTickLabel',10.^y)
hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);

%% Plot QY enhancement around gold nanorods
[xx, yy] = meshgrid(linspace( -65e-9,65e-9,2000 ),linspace( -65e-9,65e-9, 2000) );
QY = zeros(size(xx));
I_gauss = zeros(size(xx));
A_qy = scatteredInterpolant( pos(:, 1), pos(:, 2), QY_average(:)./QY_0 );

for i = 1 : numel(xx)
    QY(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_qy,ExcludedParticle );
end

figure
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],QY)
colormap jet
hold on
hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);
%% Plot near field enhancement around gold nanorods
[xx, yy] = meshgrid(linspace( -50e-9,50e-9,2000 ),linspace( -50e-9,50e-9, 2000) );
QY = zeros(size(xx));
I_gauss = zeros(size(xx));
A_exc = scatteredInterpolant( pos(:, 1), pos(:, 2), ee_full(:) );

for i = 1 : numel(xx)
    I_exc(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_exc,ExcludedParticle );
end


figure
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],I_exc)
colormap jet
hold on
hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);

%% Plot total fluorescence enhancement around gold nanorods
[xx, yy] = meshgrid(linspace( -65e-9,65e-9,2000 ),linspace( -65e-9,65e-9, 2000) );
I_fe = zeros(size(xx));
A_fe = scatteredInterpolant( pos(:, 1), pos(:, 2), ee_full(:).*QY_average(:)./0.01 );

for i = 1 : numel(xx)
    I_fe(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_fe,ExcludedParticle );
end


figure
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],I_fe)
colormap jet
hold on
hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);
%% Plot total decay rate enhancement
[row, col] = size(tot);
if col == 3
tot_avg = ( tot(:,1) + tot(:,2) + tot(:,3) ) ./3 ;
elseif col == 2
tot_avg = (tot(:,1) + 2.*tot(:, 2) )./3 ;
end 

[xx, yy] = meshgrid(linspace( -50e-9,50e-9,2000 ),linspace( -50e-9,50e-9, 2000) );
I_lifetime = zeros(size(xx));
A_lifetime = scatteredInterpolant( pos(:, 1), pos(:, 2), tot_avg(:) );

for i = 1 : numel(xx)
    I_lifetime(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_lifetime,ExcludedParticle );
end


figure
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],I_lifetime)
colormap jet
hold on
%hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);

%% Plot radiative rate enhancement 
[row, col] = size(rad);
if col == 3
rad_avg = ( rad(:,1) + rad(:,2) + rad(:,3) ) ./3 ;
elseif col == 2
rad_avg = (rad(:,1) + 2.*rad(:, 2) )./3 ;
end 

[xx, yy] = meshgrid(linspace( -50e-9,50e-9,2000 ),linspace( -50e-9,50e-9, 2000) );
I_rad = zeros(size(xx));
A_rad = scatteredInterpolant( pos(:, 1), pos(:, 2), rad_avg(:),'natural' );

for i = 1 : numel(xx)
    I_rad(i) = GetFEfromTable( xx(i),0,yy(i),pos,A_rad,ExcludedParticle );
end


figure
imagesc([min(min(xx))*1e9,-min(min(xx))*1e9],[- min(min(yy))*1e9,min(min(yy))*1e9],I_rad)
colormap jet
hold on
%hline(diameter*0.5,'w--')
p = drawNR([0,0],diameter,height);


%% Plot non-radiative rate enhancement
