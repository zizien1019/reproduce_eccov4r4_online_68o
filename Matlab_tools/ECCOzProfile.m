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

tracer_dir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';


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



%% 2. load diagnostic hb

% directory containing the file

MXLD_dir = '/Users/zengzien/Code_Research/GreatLakes/20240906_corr/MXLDEPTH_monthly/';

% Load Files

[MXLD] = rdmds([MXLD_dir 'MXLDEPTH'], [224244, 227904]);


%% 

stepsize = 3600;

figure
subplot(1,2,1)
    plot(squeeze(ptr_prof(:,1)), zc_ecco, 'r', 'LineWidth',1.5)
    hold on
    plot(squeeze(ptr_prof(:,2)), zc_ecco, 'b', 'LineWidth',1.5)
    yline(-MXLD(nx,ny,1), 'LineWidth',1.5)
    yline(-MXLD(nx,ny,2),'-.', 'LineWidth',1.5)

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
    width=500;
    height=350;
    set(gcf,'position',[x0,y0,width,height])

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');












