clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

%% 1. load tracer data

data_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';

itr_tr = 227760-730*72 : 730 : 227760-730;
ptr = rdmds([data_dir 'PTRtave01'],itr_tr);

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total

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
stress_dir = '/Users/zengzien/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_STRESS_05DEG_MONTHLY_V4R4';

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

rhodir = '20240906_corr/ziens_Rhoanoma_monthly/';
[Rho,~] = rdmds([rhodir 'RHOAnoma_mon_mean'], nan);
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

%

cor_tzave_we = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        we = squeeze(w_upward(i,j,:));
        cor_tzave_we(i,j) = corr(tzave(:), we(:));
    end
end


%% 1. load tracer data

data_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';

itr_tr = 227760-730*72 : 730 : 227760-730;
ptr = rdmds([data_dir 'PTRtave01'],itr_tr);

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total

ptr_surf = zeros(90, 1170, length(itr_tr));

for i = 1 : length(itr_tr)
    ptr_temp = ptr(:,:,1,i);
    ptr_temp(~maskc(:,:,1)) = nan;
    ptr_surf(:,:,i) = ptr_temp;
end

% clear ptr ptr_temp

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

%% 2. load diagnostic hb
% 
% % directory containing the file
% 
% MXLD_dir = '/Users/zengzien/Code_Research/GreatLakes/20240906_corr/MXLDEPTH_monthly/';
% 
% % Load Files
% 
% [MXLD, itr_hb] = rdmds([MXLD_dir 'MXLDEPTH'], nan);
% 
% itr_hb = itr_hb(240:311);
% MXLD = MXLD(:,:,240:311);

%% 2. load density to determine hb

% directory containing the file

rho_dir = '/Users/zengzien/Code_Research/GreatLakes/20240906_corr/ziens_Rhoanoma_monthly/';

% Load Files

itr = [176052 176748 177492 178212 178956 179676 180420 181164 181884 182628 183348 184092 184836 185508 186252 186972 187716 188436 189180 189924 190644 191388 192108 192852 193596 194268 195012 195732 196476 197196 197940 198684 199404 200148 200868 201612 202356 203028 203772 204492 205236 205956 206700 207444 208164 208908 209628 210372 211116 211812 212556 213276 214020 214740 215484 216228 216948 217692 218412 219156 219900 220572 221316 222036 222780 223500 224244 224988 225708 226452 227172 227903];

rho = rdmds([rho_dir 'RHOAnoma_mon_mean'], itr);


hb = zeros(90, 1170, 72);

for t = 1 : 72
    for i = 1 : 90
        for j = 1 : 1170

            r = squeeze(rho(i,j,:,t));
            drdz = gradient(r) ./ gradient(zc_ecco);
            [~,index] = max(-drdz);
            hb(i,j,t) = - zc_ecco(index);
        end
    end
end

clear rho

cor_tzave_hb = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        hb_here = squeeze(hb(i,j,:));
        cor_tzave_hb(i,j) = corr(tzave(:), hb_here);
        logcor_tzave_hb(i,j) = corr(log(tzave(:)), log(hb_here));
    end
end


%% 1. load tracer data

data_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';

itr_tr = 227760-730*72 : 730 : 227760-730;
ptr = rdmds([data_dir 'PTRtave01'],itr_tr);

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total

ptr_surf = zeros(90, 1170, length(itr_tr));

for i = 1 : length(itr_tr)
    ptr_temp = ptr(:,:,1,i);
    ptr_temp(~maskc(:,:,1)) = nan;
    ptr_surf(:,:,i) = ptr_temp;
end

% clear ptr ptr_temp

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

%% 2. load velocity for horizontal divergence

% directory containing the file

vel_dir = '/Users/zengzien/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_OCEAN_VEL_05DEG_MONTHLY_V4R4/';

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
div = zeros(90,1170);

for k = 14+1 : 14+72
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);

    u_pre = double(ncread(fullFileName, 'EVEL'));
    u = [u_pre(361:720,:,1); u_pre(1:360,:,1)];

    v_pre = double(ncread(fullFileName, 'NVEL'));
    v = [v_pre(361:720,:,1); v_pre(1:360,:,1)];

    % curl = (tau_y)_px - (tau_x)_py
    % theta = latitude !
    theta = y_podaac_t *pi/180;
    r = 6.4e6;
    
    [v_py, v_px] = gradient(v);
    [u_py, u_px] = gradient(u .* cos(theta));
    div_pd = (u_px + v_py);
    
    
    div(:,:,k-14) = interp2(x_podaac,y_podaac,div_pd.', xc_ecco,yc_ecco);

end

%%
cor_tzave_div = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        div_h = squeeze(div(i,j,:));
        cor_tzave_div(i,j) = corr(tzave(:), div_h(:));
    end
end



%%
figuredir = '/Users/zengzien/UMich Research/zien_paper_1/compound_figures_inclusive/';
tri = delaunay(xc_ecco,yc_ecco);

figure(1)
    subplot(10,4,[2 3 6 7 10 11])
    trisurf(tri,xc_ecco,yc_ecco, logcor_tzave_hb)
    % trisurf(tri,xc_ecco,yc_ecco, cor_tzave_hb)
    lighting phong
    shading interp
    colorbar EastOutside
    clim([-1 1])
    axis([0 360 -90 90])
    box on
    xlabel(['Longitude' newline newline '(a)'])
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -45 0 45 90],...
        'yticklabel',{['90' char(176) 'S'], ['45' char(176) 'S'], ['0' char(176)], ['45' char(176) 'N'], ['90' char(176) 'N']})
    colorMap = load('vik.mat');
    colormap(colorMap.vik)
    view(2)

%
    subplot(10,4,[21 22 25 26 29 30])
    trisurf(tri,xc_ecco,yc_ecco, cor_tzave_div)
    lighting phong
    shading interp
    colorbar EastOutside
    clim([-1 1])
    box on
    axis([0 360 -90 90])
    xlabel(['Longitude' newline newline '(b)'])
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -45 0 45 90],...
        'yticklabel',{['90' char(176) 'S'], ['45' char(176) 'S'], ['0' char(176)], ['45' char(176) 'N'], ['90' char(176) 'N']})
    colorMap = load('vik.mat');
    colormap(colorMap.vik)
    view(2)

%
    ax = subplot(10,4,[21 22 25 26 29 30]+2);
    trisurf(tri,xc_ecco,yc_ecco, cor_tzave_we)
    lighting phong
    shading interp
    colorbar EastOutside
    clim([-1 1])
    axis([0 360 -90 90])
    box on
    xlabel(['Longitude' newline newline '(c)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    colorMap = load('vik.mat');
    colormap(colorMap.vik)
    view(2)

    x0=10;
    y0=10;
    width=800;
    height=500;
    set(gcf,'position',[x0,y0,width,height])
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
    
    % saveas(gcf, [figuredir 'figure10_log.png'])
    % savefig(gcf, [figuredir 'figure10_log.fig'])













