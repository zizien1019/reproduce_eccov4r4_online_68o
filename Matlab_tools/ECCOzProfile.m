clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

long_interest = 230;
lat_interest = 30;

[err, wheremin] = min((xc_ecco-long_interest).^2+(yc_ecco-lat_interest).^2, [], 'all');
ny = ceil(wheremin/90);
nx = mod(wheremin, 90);

%% 1. load tracer data

% specify where tracer data are stored
tracer_dir = '';


% try take February & September

itr = [224110, 227760];
ptr = rdmds([tracer_dir 'PTRtave01'], itr);

DaysAfterStart = floor(itr/24);
DateString = datestr(727564+DaysAfterStart);

ptr_prof = zeros(50, 2);

for i = 1 : length(itr)
    ptr_temp = squeeze(ptr(nx,ny,:,i));
    ptr_temp(~maskc(nx,ny,:)) = nan;
    ptr_temp(ptr_temp<1e-7) = 1e-7;
    ptr_prof(:,i) = ptr_temp;
end

clear ptr ptr_temp

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;


%% 2. load density to determine hb

% directory containing the density anomoly files
rho_dir = '';

rho = rdmds([rho_dir 'RHOAnoma_mon_mean'], [224244, 227904]);

MLD = zeros(90, 1170, 2);

for t = 1 : 2
    for i = 1 : 90
        for j = 1 : 1170

            r = squeeze(rho(i,j,:,t));
            drdz = gradient(r) ./ gradient(zc_ecco);
            [~,index] = max(-drdz);
% definition consistent with the diagnostic MLD
% in MITgcm
            MLD(i,j,t) = zc_ecco(index);
        end
    end
end

clear rho


%% 

stepsize = 3600;

figure
    plot(squeeze(ptr_prof(:,1)), zc_ecco, 'r', 'LineWidth',1.5)
    hold on
    plot(squeeze(ptr_prof(:,2)), zc_ecco, 'b', 'LineWidth',1.5)
    yline(-MLD(nx,ny,1), 'LineWidth',1.5)
    yline(-MLD(nx,ny,2),'-.', 'LineWidth',1.5)

    ylabel('Depth (m)')
    xlabel('\tau (g/m^3)')
    % title(['Plastic concentration profile\newline' ...
    %     '\rho_p=900, d=1\mum, @30^oN, 130^oW'])
    xscale('log')
    axis([1e-4 0.2 -200 0])
    legend('\tau_{JULY}', '\tau_{DEC}', '', '', 'Location','best')
    grid on

    x0=10;
    y0=10;
    width=250;
    height=350;
    set(gcf,'position',[x0,y0,width,height])












