
clear
clc
set(0,'defaultAxesFontSize',12,'DefaultLineLineWidth',1)
set(groot,'defaultAxesXGrid','on','defaultAxesYGrid','on')

% Nstep = [219730 220460 221190 221920 222650 223380 224110 224840 225570 226300 227030 227760];
Nstep = [191260 191260];

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

% 
datadir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_2d/';
% datadir = '/Users/zengzien/Code_Research/GreatLakes/20240401_2d/ziens_2d_s/';
figuredir = '/Users/zengzien/UMich Research/zien_paper/figures/';
var_name = 'PTRtave01';
Tit_z = '';
% Tit_z = ['\rho_p=900, d=10\mum' '\newline' 'surface distribution'];
% Tit_z = ['Neutrally buoyant' '\newline' 'surface distribution'];


% datadir = '/Users/zengzien/Code_Research/GreatLakes/20240401_2d/ziens_2d_s/';
% figuredir = '/Users/zengzien/Code_Research/GreatLakes/20240401_2d/';
% var_name = 'PTRACER01';
% Tit_z = ['2D surface tracer' '\newline' 'surface distribution'];


% Tit_y = ['vertical slice at 38^oS'];
Tit_y = '';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz=1;
v_plot = ptr(:,:,nz);

% v_plot = zeros(90,1170);
% for i = 1 : 90
%     for j = 1 : 1170
%         kplot = 1;
%         for k = 1 : 49
%             if ~maskc(i,j,kplot)
%                 kplot = kplot;
%             elseif maskc(i,j,kplot) && ~maskc(i,j,kplot+1)
%                 kplot = kplot;
%             elseif maskc(i,j,kplot) && maskc(i,j,kplot+1)
%                 kplot = kplot + 1;
%             end
%         end
%         v_plot(i,j) = ptr(i,j,kplot);
%     end
% end


xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;
tri = delaunay(xc_ecco,yc_ecco);

    DaysAfterStart = floor(Nstep/24);
    DateString = datestr(727564+DaysAfterStart);

    
    figure(1)
    subplot(14,1,[2:13])
    h = trisurf(tri,xc_ecco,yc_ecco, v_plot);
    % hold on
    % tricont(xc_ecco,yc_ecco, tri, ptr(:,:,nz), [1e-1 3e-1 5e-1], 'linewidth', 2)
    axis([0 360 -90 90])

    title(Tit_z)
    
    % adjust the max and min properly
    clim([5e-4 max(v_plot,[],'all')]);
    % clim([1e-4 1])
    % clim([1e-5 10])
    
    lighting phong
    shading interp
    % colorbar northoutside
    title(colorbar,'g/m^3','FontSize',12);
    colormap("turbo")
    % colormap(slanCM(143))
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-50 0 50],...
        'yticklabel',{['50' char(176) 'S'], ['0' char(176)], ['50' char(176) 'N']})

    view(2)
    set(gca,'ColorScale','log')

    x0=10;
    y0=10;
    width=350;
    height=300;
    set(gcf,'position',[x0,y0,width,height])
    % saveas(gcf, [figuredir 'check_' var_name '_z.png'])


%% Projection onto regular grid

y_to_plot = -38;

xq = zeros(90,1170,50);
for k = 1 : 50
    xq(:,:,k) = xc_ecco;
end
yq = zeros(90,1170,50);
for k = 1 : 50
    yq(:,:,k) = yc_ecco;
end
zq = zeros(90,1170,50);
for k = 1 : 90
    zq(k,1,:) = zc_ecco(1:50);
end
for k = 2 : 1170
    zq(:,k,:) = zq(:,1,:);
end

pplease = (abs(y_to_plot-yq)<=0.5);
xt = xq(pplease);
zt = zq(pplease);
tri = delaunay(xt,zt);
    
    % interpolate
    figure(2)
    subplot(14,1,[2:13])
    hp = trisurf(tri,xt,zt,ptr(pplease));
    title(Tit_y)

    % adjust the max and min properly
    % clim([5e-4 max(ptr,[],'all')]);
    % clim([1e-6 max(ptr,[],'all')]);
    clim([1e-6 20]);
    axis([0 360 -6000 0])
    lighting phong
    shading interp
    % colorbar EastOutside
    title(colorbar,'g/m^3','FontSize',12);
    colormap("turbo")
    xlabel('Longitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ylabel('Depth (m)')
    view(2)
    set(gca,'ColorScale','log')
    % set(gca,'ColorScale','linear')

    x0=10;
    y0=10;
    width=350;
    height=200;
    set(gcf,'position',[x0,y0,width,height])
    % saveas(gcf, [figuredir 'check_' var_name '_y_6k.png'])



%%
    figure(3)
    subplot(14,1,[2:13])
    hp = trisurf(tri,xt,zt,ptr(pplease));
    title(Tit_y)

    % adjust the max and min properly
    % clim([5e-4 max(ptr,[],'all')]);
    % clim([1e-6 max(ptr,[],'all')]);
    clim([1e-6 20]);
    axis([0 360 -1000 0])
    lighting phong
    shading interp
    % colorbar EastOutside
    % title(colorbar,'g/m^3','FontSize',12);
    colormap("turbo")
    xlabel('Longitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ylabel('Depth (m)')
    view(2)
    set(gca,'ColorScale','log')
    % set(gca,'ColorScale','linear')

    x0=10;
    y0=10;
    width=350;
    height=200;
    set(gcf,'position',[x0,y0,width,height])
    % saveas(gcf, [figuredir 'check_' var_name '_y_1k.png'])
