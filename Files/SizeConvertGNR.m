%% Calculates the size and shape of GNR based on given parameters
%  Find nanorod shape parameters:
%  - of the same volume, but different aspect ratio to the given shape.
%  - of the same aspect ratio, but different volume to the given shape.
function [ h, d ] = SizeConvertGNR( height, diameter, base  )


[AR,vol] = GNRVol(height,diameter);

%% Find particles with the same volume but different AR
switch base
    case 'vol'

prompt = 'Input the aspect ratios you want to calculate ';
AR_new = input(prompt);
for i = 1 : length(AR_new)
    d(i) = ( 4/3 * vol/(AR_new(i) - 1/3))^(1/3);
end

h = d.*AR_new;
d
h
disp(['Aspect ratio = ',mat2str(AR_new)])

%% Find particles with the same AR but different volume
    case 'AR'
prompt = 'Input the volume ratios you want to calculate ';
vol_ratio = input(prompt);
for i = 1 : length(vol_ratio)
    d(i) = (vol_ratio(i))^(1/3)*diameter;
end

h = d.* AR;
d
h
disp(['volume ratio = ',mat2str(vol_ratio)])
end 

