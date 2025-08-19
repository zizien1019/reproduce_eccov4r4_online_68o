clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

%% 1. load model output

tracer_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';

itr = 227760-730*11 : 730 : 227760;
ptr = rdmds([tracer_dir 'PTRtave01'],itr);

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total

DaysAfterStart = floor(itr/24);
DateString = datestr(727564+DaysAfterStart);

ptr_surf = zeros(90, 1170, length(itr)+2);
month = [12 1 2 3 4 5 6 7 8 9 10 11 12 1];

for i = 1 : length(itr)
    ptr_temp = ptr(:,:,1,i);
    ptr_surf(:,:,i+1) = ptr_temp;
end
ptr_surf(:,:,1) = ptr_surf(:,:,13);
ptr_surf(:,:,14) = ptr_surf(:,:,2);

% concentration to number density
% 1 #/km^2
particle_mass = 900 * 4/3*pi*(10e-6/2)^3; % kg/particle
ptr_surf = 1e3 * ptr_surf .* hFacC(:,:,1) * DRF(1) / particle_mass; % #/km^2

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

[~,month_peak] = max(ptr_surf(:,:,2:13),[],3,"includemissing");
peak_time_ecco = month(month_peak);
peak_time_ecco(~maskc(:,:,1)) = nan;
peak_time_ecco(abs(yc_ecco)>35) = nan;

save('peak_time_ecco');

%% 2. load CYGNSS data

plast_dir = '/Users/zengzien/Downloads/ECCO_V4r4_PODAAC/CYGNSS_L3_MICROPLASTIC_V1.0/';

% Load Files
% Specify the folder where the files live.
myFolder = plast_dir;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.nc'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

baseFileName = theFiles(1).name;
fullFileName = fullfile(theFiles(1).folder, baseFileName);
lon_cygnss = double(ncread(fullFileName, 'lon'));
lat_cygnss = double(ncread(fullFileName, 'lat'));
[x_cygnss,y_cygnss] = ndgrid(lon_cygnss,lat_cygnss);

ptr_temp = zeros(1440, 297);
ptr_cygnss = zeros(90, 1170, 14);
rising_time_cygnss = zeros(90, 1170);

month_cygnss = [4 5 6 7 8 9 10 11 12 1 2 3];

for k = 1 : 12
    MP = double(ncread(fullFileName, 'MP_concentration'));
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    ptr_temp = mean(MP,3,"includemissing");
    for i = 1 : 90
        for j = 1 : 1170
            [~, nx] = min(abs(xc_ecco(i,j)-lon_cygnss));
            [~, ny] = min(abs(yc_ecco(i,j)-lat_cygnss));

            ptr_cygnss(i,j,month_cygnss(k)+1) = mean( ptr_temp(...
                max(nx-2,1):min(nx+2,1440), ...
                max(ny-2,1):min(ny+2,297) ...
                ), 'all', 'omitmissing');
        end
    end
end
ptr_cygnss(:,:,1) = ptr_cygnss(:,:,13);
ptr_cygnss(:,:,14) = ptr_cygnss(:,:,2);
month = [12 1 2 3 4 5 6 7 8 9 10 11 12 1];

[~,month_peak] = max(ptr_cygnss(:,:,2:13),[],3,"includemissing");
peak_time_cygnss = month(month_peak);
peak_time_cygnss(~maskc(:,:,1)) = nan;
peak_time_cygnss(abs(yc_ecco)>35) = nan;

save('peak_time_cygnss');

%% 3. calculate delta t_p

clear
clc

load('XYZ_ecco.mat');
xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;
tri = delaunay(xc_ecco,yc_ecco);
load('Geo_ecco.mat');
load("peak_time_cygnss.mat");
load("peak_time_ecco.mat")

% calculate \Delta T_p

del_t = peak_time_cygnss - peak_time_ecco ;
del_t(del_t < -6) = del_t(del_t < -6) + 12;
del_t(del_t >  6) = del_t(del_t >  6) - 12;

subplot(3,1,2)
    h = histogram(del_t(:), -6.5:6.5);
    percentage = sum( h.BinCounts(1+3:end-3) ) / sum(h.BinCounts);
    axis([-6.5 6.5 -inf inf])
    xlabel('$\Delta T_p$ $\mbox{(c)}$', 'Interpreter','latex', 'FontName','Times New Roman')
    xticks(-6:6)
    ylabel('# of points')
    mu = num2str( mean(del_t, 'all', 'omitmissing') );
    si = num2str( std(del_t(:), 'omitmissing') );
    text(-5.5,2100,['Mean = ' mu newline 'Std = ' si])

    view(2)
    x0=600;
    y0=50;
    width=600;
    height=500;
    set(gcf,'position',[x0,y0,width,height])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');














