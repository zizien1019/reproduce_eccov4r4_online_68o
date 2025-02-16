
% Locate CYGNSS data directory
cyg1 = '/Users/zengzien/Downloads/ECCO_V4r4_PODAAC/';
cyg2 = 'CYGNSS_L3_MICROPLASTIC_V1.0/';
data_dir = [cyg1 cyg2];

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(data_dir, '*.nc'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% load the files
% ncdisp
% ncread Syntax:
% vardata = ncread(source,varname)
% vardata = ncread(source,varname,start,count)
% vardata = ncread(source,varname,start,count,stride)

baseFileName = theFiles(1).name;
fullFileName = fullfile(theFiles(1).folder, baseFileName);
lat = double(ncread(fullFileName, 'lat'));
lon = double(ncread(fullFileName, 'lon'));

% for k = 1 : Timespan
for k = 1
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    MiPla{k,1} = double(ncread(fullFileName, 'MP_concentration'));

end

% cyg.ddmi.s20170402-000000-e20170430-000000.l3.grid-microplastic.a10.d10.nc;

surf(lon,37:-0.25:-37,log10(MiPla{1}(:,:,1)).','EdgeColor','none')
title('CYGNSS Ocean Microplastics Data at 2017-Apr-02, in log(#)/km^2')
lighting phong
shading interp
colorbar EastOutside
colormap("jet")
xlabel('Longitude(^o)')
ylabel('Latitude(^o)')
view(2)

