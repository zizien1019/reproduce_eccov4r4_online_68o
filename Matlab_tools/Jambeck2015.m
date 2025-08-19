clear
clc
set(0,'defaultAxesFontSize',12,'DefaultLineLineWidth',1)
set(groot,'defaultAxesXGrid','on','defaultAxesYGrid','on')

% Load data
% Jambeck 2015 :
% total mismanaged waste ~ 3.1865e7 tons in 2010
load("borders_dataset.mat");
load("waste_gdp.mat");
waste = waste_gdp(:,1);
gdp = waste_gdp(:,2);
load('XYZ_ecco.mat');
load('Geo_ecco.mat')
rA = RAC;
dz = ones(90,1170) * DRF(1);

waste_adjusted = waste;

% Projection onto ecco grid
% find the closest 1:2 values
% for an evenly distributed grid points
% 90 * 1170 = 150 * 351
xb_0 = linspace(-180,180,360);
yb_0 = linspace(-90,90,180);
[xq_0,yq_0] = ndgrid(xb_0,yb_0);
vq_0 = zeros(size(xq_0));
xq_0 = xq_0(:); yq_0 = yq_0(:); vq_0 = vq_0(:);

vc = zeros(size(xc_ecco));
n_extra = zeros(302);

% 225: Indonesia
% 21: Brazil 
% 24: Canada
% 209: United States
% 175: Russia
for i = 1 : 302
% acquire extra query points country by country
    x_ext = lon{i}.';
    y_ext = lat{i}.';


    if i == 24
        x_ext(y_ext > 54) = [];
        y_ext(y_ext > 54) = [];
    elseif i == 209
        x_ext(y_ext > 60) = [];
        y_ext(y_ext > 60) = [];
    elseif i == 175
        x_ext(y_ext > 70) = [];
        y_ext(y_ext > 70) = [];
    end


    v_ext = ones(size(x_ext)) * waste_adjusted(i);
% add data to query points country by country
    xq = [xq_0; x_ext; x_ext+0.5; x_ext+0.5; x_ext-0.5; x_ext-0.5; x_ext+1; x_ext+1; x_ext-1; x_ext-1];
    yq = [yq_0; y_ext; y_ext+0.5; y_ext-0.5; y_ext+0.5; y_ext-0.5; y_ext+1; y_ext-1; y_ext+1; y_ext-1];
    vq = [vq_0; v_ext; v_ext; v_ext; v_ext; v_ext; v_ext; v_ext; v_ext; v_ext];
% clean the nan's
    vq(isnan(xq))=[];
    yq(isnan(xq))=[];
    xq(isnan(xq))=[];
% interpolate
    F = scatteredInterpolant(xq,yq,vq,'natural','nearest');
    vc_inter = F(xc_ecco,yc_ecco);
    vc_inter(vc_inter<0) = 0;
% mask out the land
    vc_inter(~maskc(:,:,1)) = 0;
% normalize the amount of tracers
% evenly distributed along the coast of each country
    pop_up = (vc_inter~=0);
    n_extra = sum(sum(double(pop_up)));
    if n_extra~=0
        vc_inter = vc_inter/n_extra;
    end
% update for next iteration
    vc = vc + vc_inter;
end

vc = vc * 3.1865e7 / sum(vc,'all');
vc = vc * 0.4;




%% Convert to source input format

vc(isnan(vc))=0;
load('Geo_ecco.mat')
dz = hFacC(:,:,1) * DRF(1);
Stons = vc; % tons
Sppmps = vc ./ RAC ./ dz  * 1e6 / 365 / 86400; % ppm per s




%% Visualize

clear
clc
set(0,'defaultAxesFontSize',16,'DefaultLineLineWidth',1)
set(groot,'defaultAxesXGrid','on','defaultAxesYGrid','on')

% Load data
% Jambeck 2015 :
% total mismanaged waste ~ 3.1865e7 tons in 2010
load("borders_dataset.mat");
load("waste_gdp.mat");
waste = waste_gdp(:,1);
gdp = waste_gdp(:,2);
load('XYZ_ecco.mat');
load('Geo_ecco.mat')
rA = RAC;
dz = ones(90,1170) * DRF(1);

load("Q_ppmps.mat");
sf = Sppmps(:,:,1);
sf(~maskc(:,:,1)) = nan;
sf(sf==0)=nan;

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;
tri = delaunay(xc_ecco,yc_ecco);

figure(2)
    subplot(14,1,2:13)
    h = trisurf(tri,xc_ecco,yc_ecco,sf);
    axis([0 360 -90 90])
    clim([1e-11 3e-8])
    lighting phong
    shading interp
    colorbar EastOutside
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -45 0 45 90],...
        'yticklabel',{['90' char(176) 'S'], ['45' char(176) 'S'], ['0' char(176)], ['45' char(176) 'N'], ['90' char(176) 'N']})
    box on
    title(colorbar,'g m^{-3} s^{-1}');
    set(gca,'ColorScale','log')
    view(2)

    x0=10;
    y0=10;
    width=500;
    height=300;
    set(gcf,'position',[x0,y0,width,height])
    movegui(gca,'center')

    % figuredir = '/Users/zengzien/UMich Research/zien_paper/figures/';
    figuredir = './';
    saveas(gcf, [figuredir 'waste_40.png'])
    savefig([figuredir 'waste_40.fig'])















































