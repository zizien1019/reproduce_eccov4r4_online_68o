clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

%% 1. load tracer data

data_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr/';

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total
itr_tr = 227760-730*71 : 730 : 227760;
ptr = rdmds([data_dir 'PTRtave01'],itr_tr);

ptr_surf = zeros(90, 1170, length(itr_tr));

for i = 1 : length(itr_tr)
    ptr_temp = ptr(:,:,1,i);
    ptr_temp(~maskc(:,:,1)) = nan;
    ptr_surf(:,:,i) = ptr_temp;
end

% clear ptr ptr_temp

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

%% 2. load wind stress for upwelling

% directory containing the file
% the Ekman velocity is computed from the surface wind stress data
% also available on PODAAC
stress_dir = '~/ECCO_V4r4_PODAAC/ECCO_L4_STRESS_05DEG_MONTHLY_V4R4';

% Load Files
% Specify the folder where the files live.
myFolder = stress_dir;

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
stressFiles = dir(filePattern);

baseFileName = stressFiles(1).name;
    fullFileName = fullfile(stressFiles(1).folder, baseFileName);
    x_podaac = double(ncread(fullFileName, 'longitude'));
    x_podaac_t = [x_podaac(361:720,1); x_podaac(1:360,1)+360];
    y_podaac_t = double(ncread(fullFileName, 'latitude'));
    y_podaac_t = y_podaac_t.';
    [x_podaac, y_podaac] = meshgrid(x_podaac_t(:), y_podaac_t(:));

rho0 = 1025;
w_upwelling_pd = zeros(720,360);
w_upward = zeros(size(ptr_surf));


for k = 240+1 : 240+72
    baseFileName = stressFiles(k).name;
    fullFileName = fullfile(stressFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);

    tau_x_pre = double(ncread(fullFileName, 'oceTAUE'));
    tau_x = [tau_x_pre(361:720,:); tau_x_pre(1:360,:)];

    tau_y_pre = double(ncread(fullFileName, 'oceTAUN'));
    tau_y = [tau_y_pre(361:720,:); tau_y_pre(1:360,:)];

    % curl = (tau_y)_px - (tau_x)_py
    % theta = latitude !
    theta = y_podaac_t *pi/180;
    r = 6.4e6;
    
    [~, tau_y_px] = gradient(tau_y);
    [tau_x_py, ~] = gradient(tau_x .* cos(theta));
    curl = (tau_y_px - tau_x_py) / (0.5 *pi/180) ./ (r*cos(theta));
    
    f = 2 * 7.2921e-5 * sin(theta);
    
    w_upwelling_pd = curl ./ f / rho0;
    % F = scatteredInterpolant(x_podaac(:),y_podaac(:),w_upwelling_pd(:),'natural');
    % w_upwelling(:,:,k-240) = F(xc_ecco,yc_ecco);
    w_upward(:,:,k-240) = interp2(x_podaac,y_podaac,w_upwelling_pd.', xc_ecco,yc_ecco);

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

cor_tzave_we = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        we = squeeze(w_upward(i,j,:));
        cor_tzave_we(i,j) = corr(tzave(:), we(:));
    end
end


%%
tri = delaunay(xc_ecco,yc_ecco);
% specify where to store the figures
figuredir = '';
figure(1)
    subplot(14,1,[2:13])
trisurf(tri,xc_ecco,yc_ecco, cor_tzave_we)
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
    saveas(gcf, [figuredir '10mum_corr_we_tsurf.png'])