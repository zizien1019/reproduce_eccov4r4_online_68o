%% interpolate >> for source from other grid
RAC = rdmds('RAC');
rA = RAC;
DRF = rdmds('DRF');
dz = ones(90,1170) * DRF(1);
load('XYZ_ecco.mat');
Source = zeros(90,1170);
variable_to_write = zeros(90,1170,50);

% get source information
lat = plasticWaste(2:end, 1);
lon = plasticWaste(2:end, 2);
src = plasticWaste(2:end, 3);

% interpolate to generic grid
x = xc_ecco;
x(~maskc(:,:,1)) = 200;
y = yc_ecco;
y(~maskc(:,:,1)) = 100;
for i = 1 : length(lon)
% for i = 1 : 1
    distance = (x-lon(i)).^2 + (y-lat(i)).^2;
    [M,I] = min(distance,[],'all');
%     Source(I) = src(i);
    Source(I) = src(i) / rA(I) / dz(I) ; % tons >> ton/m^3
end
Source = Source * 10^6 ; % ton/m^3 = 10^3 kg/m^3 = 10^6 g/m^3
Source = Source / 365 / 86400 ; % adjust to every second

% write binary file
variable_to_write(:,:,1) = Source;
fileID = fopen('myfile.bin','w');
fwrite(fileID,variable_to_write,'single','ieee-be');
fclose(fileID);

%% generate source >> based on generic grid
RAC = rdmds('ecco_geometry/RAC');
rA = RAC;
DRF = rdmds('ecco_geometry/DRF');
dz = ones(90,1170) * DRF(1);
load('XYZ_ecco.mat');
Source = zeros(90,1170);
variable_to_write = zeros(90,1170,50);

%  source
Source = (yc_ecco<40) & (yc_ecco>-40);
Source = single(Source);
Source = Source * 1e-3;

%% generate source >> for advection_in_gyre test
% in 60 * 60 * 10 domain
% variable_to_write = zeros(90,40,15);
variable_to_write = 1e-4 * ones(90,40,15);
% Source = ones(90,40)*1e-4;

% write binary file
% variable_to_write(:,:,1) = Source;
fileID = fopen('surf1_adv_test.bin','w');
fwrite(fileID,variable_to_write,'single','ieee-be');
fclose(fileID);

%% Visualization
index = (Source~=0);
geoscatter(y(index),x(index),Source(index)/30000,Source(index),'filled');
colormap('copper')
colorbar
legend('plastic release in 2010 (tons)')
title(['Source distribution' '\newline' 'based on Jambeck et al. 2015'])







