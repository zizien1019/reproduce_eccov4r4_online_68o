%% 1. load averaged tracer data

clear
clc

load('ecco_geometry/Geo_ecco.mat');
load('XYZ_ecco.mat');
tri = delaunay(xc_ecco,yc_ecco);

% specify where the data are stored
datadir = ''; 

% iteration numbers corresbonding to 2012 Jan - 2017 Dec
itr = 227760-730*72 : 730 : 227760-730;

% read 6-years of monthly-averaged surface tracer value
% so 72 records in total
ptr_surf = zeros(90, 1170, 72);
% in our tests PTRACER01 is the label for the particles
% with density 900 kg m-3 and size 10 \mu m
[ptr] = rdmds([datadir 'PTRtave01'], itr);

Timespan = 1 : 72;
for i = Timespan

    ptr_temp = ptr(:,:,1,i);
    ptr_temp(~maskc(:,:,1)) = nan;
    ptr_surf(:,:,i) = ptr_temp;
end

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

%% 2. load velocity for horizontal divergence

% directory containing the velocity files
% take monthly averaged ECCO velocity data
vel_dir = '~/ECCO_V4r4_PODAAC/ECCO_L4_OCEAN_VEL_05DEG_MONTHLY_V4R4/';

% Load Files

% Specify the folder where the files live.
myFolder = vel_dir;

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
    x_podaac = double(ncread(fullFileName, 'longitude'));
    x_podaac_t = [x_podaac(361:720,1); x_podaac(1:360,1)+360];
    y_podaac_t = double(ncread(fullFileName, 'latitude'));
    y_podaac_t = y_podaac_t.';
    [x_podaac, y_podaac] = meshgrid(x_podaac_t(:), y_podaac_t(:));

div_pd = zeros(720,360);
div = zeros(size(ptr_surf));

% k index are flexible
% make sure the files loaded correspond to the same time period
% (the same month) taken for the tracer values
for k = 14+1 : 14+72
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);

    u_pre = double(ncread(fullFileName, 'EVEL'));
    u = [u_pre(361:720,:,1); u_pre(1:360,:,1)];

    v_pre = double(ncread(fullFileName, 'NVEL'));
    v = [v_pre(361:720,:,1); v_pre(1:360,:,1)];

    % curl = (tau_y)_px - (tau_x)_py
    % theta = latitude
    theta = y_podaac_t *pi/180;
    r = 6.4e6;
    % the gradient operator in spherical coordinate
    [v_py, v_px] = gradient(v);
    [u_py, u_px] = gradient(u .* cos(theta));
    div_pd = (u_px + v_py);
    
    
    div(:,:,k-14) = interp2(x_podaac,y_podaac,div_pd.', xc_ecco,yc_ecco);

end


%% plot to check

% ntime = 2;
% % trisurf(tri,xc_ecco,yc_ecco, squeeze(w_upwelling(:,:,ntime)))
% % trisurf(tri,xc_ecco,yc_ecco, squeeze(ptr_zave(:,:,ntime)))
% surf(x_podaac,y_podaac,w_upwelling_pd.')
% lighting phong
% shading interp
% colorbar EastOutside
% clim([-1e-5 1e-5])
% colormap(slanCM(104))
% view(2)

%%
% 3. compute correlation
% along the third dimension (time)!
% syntax : correlation = corr([1; 2; 3], [30; 90; 80]);
% input matrices must be in columns

cor_tzave_div = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        div_h = squeeze(div(i,j,:));
        cor_tzave_div(i,j) = corr(tzave(:), div_h(:));
    end
end


%%
tri = delaunay(xc_ecco,yc_ecco);
% specify where to store the figures
figuredir = '';
figure(1)
    subplot(14,1,[2:13])
trisurf(tri,xc_ecco,yc_ecco, cor_tzave_div)
    lighting phong
    shading interp
    colorbar EastOutside
    clim([-1 1])
    axis([0 360 -90 90])
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-50 0 50],...
        'yticklabel',{['50' char(176) 'S'], ['0' char(176)], ['50' char(176) 'N']})
    % 
    colormap(slanCM(104))
    view(2)
    x0=10;
    y0=10;
    width=350;
    height=200;
    set(gcf,'position',[x0,y0,width,height])
    saveas(gcf, [figuredir '10mum_corr_div_tsurf.png'])