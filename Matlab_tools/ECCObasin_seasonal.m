clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');
xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

long_interest = 230;
lat_interest = 30;
    
[err, wheremin] = min( ((xc_ecco-long_interest).^2+(yc_ecco-lat_interest).^2) , [], 'all');
ny = ceil(wheremin/90);
nx = mod(wheremin, 90);

%% 1. load tracer data

data_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';

itr = 227760-730*23 : 730 : 227760;
ptr = rdmds([data_dir 'PTRtave01'],itr);

% monthly average
% 2012 Jan - 2017 Dec
% 72 records in total

DaysAfterStart = floor(itr/24);
DateString = datestr(727564+DaysAfterStart);

ptr_surf =  squeeze(ptr(nx,ny,1,:));

% clear ptr

ptr_ML = zeros(24,1);


%% 2. load density to determine hb

% directory containing the file

rho_dir = '/Users/zengzien/Code_Research/GreatLakes/20240906_corr/ziens_Rhoanoma_monthly/';

% Load Files

itr = [ ...
    210372;
    211116;
    211812;
    212556;
    213276;
    214020;
    214740;
    215484;
    216228;
    216948;
    217692;
    218412;
    219156;
    219900;
    220572;
    221316;
    222036;
    222780;
    223500;
    224244;
    224988;
    225708;
    226452;
    227172
    ];

itr = itr.';
rho = rdmds([rho_dir 'RHOAnoma_mon_mean'], itr);


hb = zeros(90, 1170, 24);

for t = 1 : 24
    for i = 1 : 90
        for j = 1 : 1170

            r = squeeze(rho(i,j,:,t));
            drdz = gradient(r) ./ gradient(zc_ecco);
            [~,index] = max(-drdz);
            hb(i,j,t) = zc_ecco(index);
            % ptr_ML(t) = squeeze(ptr(nx,ny,1:index)).' * squeeze(DRF(1:index));
        end
    end
end

clear rho

%% 2. load diagnostic hb

% % directory containing the file
% 
% MXLD_dir = '/Users/zengzien/Code_Research/GreatLakes/20240906_corr/MXLDEPTH_monthly/';
% 
% % Load Files
% 
% [~, itr] = rdmds([MXLD_dir 'MXLDEPTH'], nan);
% MXLD = rdmds([MXLD_dir 'MXLDEPTH'], itr(312-23:312));
% 
% ptr_ML = zeros(24,1);
% 
% for t = 1 : 24
%     dh = zeros(50,1);
%     H = MXLD(nx,ny,t);
%     for k = 1 : 50
%         if H > 0
%             dh(k) = min(H, DRF(k));
%             H = H - DRF(k);
%         elseif (H <= 0)
%             dh(k) = 0;
%         end
%     end
% 
%     ptr_ML(t) = squeeze(ptr(nx,ny,:,t)).' * dh;
% end



%% 3. load velocity for horizontal divergence

% directory containing the file

vel_dir = '/Users/zengzien/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_OCEAN_VEL_05DEG_MONTHLY_V4R4/';

% Load Files
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
div = zeros(90,1170,24);

for k = 14+72-23 : 14+72
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

    [v_py, v_px] = gradient(v ./ cos(theta));
    [u_py, u_px] = gradient(u);
    div_pd = (u_px + v_py) / r;


    div(:,:,k-14-72+24) = interp2(x_podaac,y_podaac,div_pd.', xc_ecco,yc_ecco);

end

%% 4. load wind stress for upwelling

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

w_upwelling_pd = zeros(720,360);
w_upwelling = zeros(90, 1170, 24);
rho0 = 1025;


for k = 240+72-23 : 240+72
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
    % w_upwelling(:,:,k-240-72+24) = w_upwelling_pd;

    w_upwelling(:,:,k-240-72+24) = interp2(x_podaac,y_podaac,w_upwelling_pd.', xc_ecco,yc_ecco);
    

    % F = scatteredInterpolant(x_podaac(:),y_podaac(:),w_upwelling_pd(:),'natural');
    % w_upwelling(:,:,k-240) = F(xc_ecco,yc_ecco);
    % w_upward(:,:,k-240) = interp2(x_podaac,y_podaac,w_upwelling_pd.', xc_ecco,yc_ecco);

end

%% construct mask for each basin

% Npac = ((yc_ecco>=18)&(yc_ecco<=38)&(xc_ecco>=150)&(xc_ecco<=230)) ...
%     & maskc;
% 
% ExNpac = ((yc_ecco>=5)&(yc_ecco<=60)&(xc_ecco>=150)&(xc_ecco<=230)) ...
%     & maskc;
% 
% Spac = ((yc_ecco>=-38)&(yc_ecco<=-18)&(xc_ecco>=180)&(xc_ecco<=260)) ...
%     & maskc;
% 
% Natl = ((yc_ecco>=0)&(yc_ecco<=67)&(xc_ecco>=260)&(xc_ecco<=360)) ...
%     & (~((yc_ecco>=0)&(yc_ecco<=15.5)&(xc_ecco>=260)&(xc_ecco<=278))) ...
%     & maskc;
% 
% Satl = (yc_ecco<=0) & ((xc_ecco>=278)|(xc_ecco<=20.5)) ...
%     & maskc;
% 
% Indi = ((xc_ecco>=20.5)&(yc_ecco<=26)&(xc_ecco<=140.5)) ...
%     & (~((yc_ecco>=-6.5)&(yc_ecco<=26)&(xc_ecco>=103.5)&(xc_ecco<=140.5))) ...
%     & maskc;
% 
% Medi = (yc_ecco>=30)&(yc_ecco<=50)&(xc_ecco>=0)&(xc_ecco<=50) ...
%     & maskc;


%% hor_ave

Nstps = 1 : 24;
Monthnum = 1 : 12;

figuredir = '/Users/zengzien/UMich Research/zien_paper_1/compound_figures_inclusive/';

figure(1)
    
    ax = subplot(2,2,1);
    ax.Position(1) = ax.Position(1);
    ax.Position(2) = ax.Position(2) + 0.08;
    ax.Position(4) = ax.Position(4) * 0.6;
    plot(Nstps, squeeze( - hb(nx,ny,Nstps)), 'LineWidth',1.5)
    set(gca,'xtick',1:24,...
        'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'})
    xlabel(['Month' newline newline '(a)'])
    ylabel('MLD (m)')
    grid on
    

    ax = subplot(2,2,2);
    ax.Position(1) = ax.Position(1) - 0.03;
    ax.Position(2) = ax.Position(2) + 0.08;
    ax.Position(4) = ax.Position(4) * 0.6;
    plot(Nstps, squeeze(div(nx,ny,Nstps)), 'LineWidth',1.5)
    set(gca,'xtick',1:24,...
        'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'})
    xlabel(['Month' newline newline '(b)'])
    ylabel('\nabla_h \cdot u_h (1/s)')
    grid on
    
    
    ax = subplot(2,2,3);
    ax.Position(1) = ax.Position(1);
    ax.Position(2) = ax.Position(2) + 0.15;
    ax.Position(4) = ax.Position(4) * 0.6;
    plot(Nstps, squeeze(w_upwelling(nx,ny,Nstps)), 'LineWidth',1.5)
    set(gca,'xtick',1:24,...
        'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'})
    xlabel(['Month' newline newline '(c)'])
    ylabel('w_E (m/s)')
    grid on

    ax = subplot(2,2,4);
    ax.Position(1) = ax.Position(1) - 0.03;
    ax.Position(2) = ax.Position(2) + 0.15;
    ax.Position(4) = ax.Position(4) * 0.6;
    ptr_surf(ptr_surf==0) = nan;
    plot(Nstps, squeeze(ptr_surf(Nstps)), 'LineWidth',1.5)
    set(gca,'xtick',1:24,...
        'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'})
    xlabel(['Month' newline newline '(d)'])
    ylabel('\tau (g/m^3)')
    grid on

    x0=10;
    y0=10;
    width=1000;
    height=600;
    set(gcf,'position',[x0,y0,width,height])



    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

    % saveas(gcf, [figuredir 'figure9.png'])
    % savefig(gcf, [figuredir 'figure9.fig'])









