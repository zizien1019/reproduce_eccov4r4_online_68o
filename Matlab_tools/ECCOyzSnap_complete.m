%%

% clear
% clc
% 
% color_ends_z = [1e-5 5];
% color_ends_y = [1e-6 1];
% 
% Nstep = [183230 183960 184690 185420 186150 186880 187610 188340 189070 189800 190530 191260];
% 
% load('XYZ_ecco.mat');
% load('Geo_ecco.mat');
% % Mesh used for surface plot
%     xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;
%     tri_z = delaunay(xc_ecco,yc_ecco);
% % Mesh used for slice plot
%     y_to_plot = -38;
%     xq = zeros(90,1170,50);
%     for k = 1 : 50
%         xq(:,:,k) = xc_ecco;
%     end
%     yq = zeros(90,1170,50);
%     for k = 1 : 50
%         yq(:,:,k) = yc_ecco;
%     end
%     zq = zeros(90,1170,50);
%     for k = 1 : 90
%         zq(k,1,:) = zc_ecco(1:50);
%     end
%     for k = 2 : 1170
%         zq(:,k,:) = zq(:,1,:);
%     end
%     pplease = (abs(y_to_plot-yq)<=0.5);
%     xt = xq(pplease);
%     zt = zq(pplease);
%     tri_y = delaunay(xt,zt);
% 
% % 
% datadir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_2d/';
% figuredir = '/Users/zengzien/UMich Research/zien_paper_1/compound_figures_inclusive/';
% Tit_z = '';
% Tit_y = '';



%% Figure 2

% figure(200)
% 
% var_name = 'PTRtave01';
% 
% ptr_all = rdmds([datadir var_name], Nstep);
% ptr = mean(ptr_all,4,"omitmissing");
% ptr(~maskc) = nan;
% 
% nz = 1;
% v_plot = ptr(:,:,nz);
% 
% %    
%     ax = subplot(13,2,[1 3 5 7 9 11]);
%     h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
%     axis([0 360 -90 90])
%     clim(color_ends_z)
%     title(Tit_z)
%     lighting phong
%     shading interp
%     colorMap = load('batlow.mat');
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline '(a)'])
%     ylabel('Latitude')
%     yline(-38,'-.w','LineWidth',1)
%     box on
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     set(gca,'ytick',[-90 -38 0 38 90],...
%         'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
%     view(2)
%     set(gca,'ColorScale','log')
% 
% % 
%     ax = subplot(13,2,[17 19 21 23]);
%     hp = trisurf(tri_y,xt,zt,ptr(pplease));
%     title(Tit_y)
%     clim(color_ends_y);
%     axis([0 360 -1000 0])
%     lighting phong
%     shading interp
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline '(b)'])
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     ylabel('Depth (m)')
%     view(2)
%     set(gca,'ColorScale','log')
%     box on
% 
% 
% var_name = 'PTRtave01';
% 
% ptr_all = rdmds([datadir var_name], Nstep);
% ptr = mean(ptr_all,4,"omitmissing");
% ptr(~maskc) = nan;
% 
% nz = 1;
% v_plot = ptr(:,:,nz);
% 
% %
%     ax = subplot(13,2,[1 3 5 7 9 11]+1);
%     h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
%     axis([0 360 -90 90])
%     clim(color_ends_z)
%     title(Tit_z)
%     lighting phong
%     shading interp
%     colorMap = load('batlow.mat');
%     colormap(colorMap.batlow)
%     xlabel(['Longitude'])
%     ylabel('Latitude')
%     % yline(-38,'-.w','LineWidth',1)
%     box on
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     set(gca,'ytick',[-90 -38 0 38 90],...
%         'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
%     view(2)
%     set(gca,'ColorScale','log')
% 
%     % Create shared colorbar for first row (ax1 and ax2)
%     cb1 = colorbar('Position', [0.92, 0.558, 0.02, 0.368]); % [left, bottom, width, height]
%     title(cb1, 'g m^{-3}')
% 
% % 
%     ax = subplot(13,2,[17 19 21 23]+1);
%     hp = trisurf(tri_y,xt,zt,ptr(pplease));
%     title(Tit_y)
%     clim(color_ends_y);
%     axis([0 360 -1000 0])
%     lighting phong
%     shading interp
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline '(d)'])
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     ylabel('Depth (m)')
%     view(2)
%     set(gca,'ColorScale','log')
% 
%     x0  =   10;
%     y0  =   10;
%     width   =   1300;
%     height  =   800;
%     set(gcf,'position',[x0,y0,width,height])   
% 
% 
%     % Create shared colorbar for second row (ax3 and ax4)
%     cb2 = colorbar('Position', [0.92, 0.172, 0.02, 0.24]);
%     title(cb2, 'g m^{-3}')
% 
% 
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
%     set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
% 
%     saveas(gcf, [figuredir 'figure2.png'])
%     savefig(gcf, [figuredir 'figure2.fig'])

%%

clear
clc



color_ends_z = [9e-6 1.1];
color_ends_y = [1e-6 1];

Nstep = [219730 220460 221190 221920 222650 223380 224110 224840 225570 226300 227030 227760];

load('XYZ_ecco.mat');
load('Geo_ecco.mat');
% Mesh used for surface plot
    xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;
    tri_z = delaunay(xc_ecco,yc_ecco);
% Mesh used for slice plot
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
    tri_y = delaunay(xt,zt);

% 
datadir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';
figuredir = '/Users/zengzien/UMich Research/zien_paper_1/compound_figures_inclusive/';
Tit_z = '';
Tit_y = '';



%% Figure 3
%
% figure(300)
% 
% var_name = 'PTRtave05';
% ptr_all = rdmds([datadir var_name], Nstep);
% ptr = mean(ptr_all,4,"omitmissing");
% ptr(~maskc) = nan;
% 
% nz = 1;
% v_plot = ptr(:,:,nz);
% 
% %    
%     ax = subplot(13,2,[1 3 5 7 9 11]);
%     h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
%     axis([0 360 -90 90])
%     clim(color_ends_z)
%     title(Tit_z)
%     lighting phong
%     shading interp
%     colorMap = load('batlow.mat');
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline newline '(a)'])
%     ylabel('Latitude')
%     yline(-38,'-.w','LineWidth',1)
%     box on
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     set(gca,'ytick',[-90 -38 0 38 90],...
%         'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
%     view(2)
%     set(gca,'ColorScale','log')
% 
% % 
%     ax = subplot(13,2,[17 19 21 23]);
%     hp = trisurf(tri_y,xt,zt,ptr(pplease));
%     title(Tit_y)
%     clim(color_ends_y);
%     axis([0 360 -1000 0])
%     lighting phong
%     shading interp
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline newline '(b)'])
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     ylabel('Depth (m)')
%     view(2)
%     set(gca,'ColorScale','log')
%     box on
% 
% 
% var_name = 'PTRtave05';
% 
% ptr_all = rdmds([datadir var_name], Nstep);
% ptr = mean(ptr_all,4,"omitmissing");
% ptr(~maskc) = nan;
% 
% nz = 1;
% v_plot = ptr(:,:,nz);
% 
% %
%     ax = subplot(13,2,[1 3 5 7 9 11]+1);
%     h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
%     axis([0 360 -90 90])
%     clim(color_ends_z)
%     title(Tit_z)
%     lighting phong
%     shading interp
%     colorMap = load('batlow.mat');
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline newline '(a)'])
%     ylabel('Latitude')
%     yline(-38,'-.w','LineWidth',1)
%     box on
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     set(gca,'ytick',[-90 -38 0 38 90],...
%         'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
%     view(2)
%     set(gca,'ColorScale','log')
% 
%     % Create shared colorbar for first row (ax1 and ax2)
%     cb1 = colorbar('Position', [0.92, 0.558, 0.02, 0.368]); % [left, bottom, width, height]
%     title(cb1, 'g m^{-3}')
% 
% % 
%     ax = subplot(13,2,[17 19 21 23]+1);
%     hp = trisurf(tri_y,xt,zt,ptr(pplease));
%     title(Tit_y)
%     clim(color_ends_y);
%     axis([0 360 -1000 0])
%     lighting phong
%     shading interp
%     colormap(colorMap.batlow)
%     xlabel(['Longitude' newline newline '(b)'])
%     set(gca,'xtick',[0 90 180 270 360],...
%         'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
%     ylabel('Depth (m)')
%     view(2)
%     set(gca,'ColorScale','log')
% 
%     x0  =   10;
%     y0  =   10;
%     width   =   1300;
%     height  =   800;
%     set(gcf,'position',[x0,y0,width,height])   
% 
% 
%     % Create shared colorbar for second row (ax3 and ax4)
%     cb2 = colorbar('Position', [0.92, 0.172, 0.02, 0.24]);
%     title(cb2, 'g m^{-3}')
% 
% 
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
%     set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
% 
%     saveas(gcf, [figuredir 'figure3.png'])
%     savefig(gcf, [figuredir 'figure3.fig'])

%% Figure 4

figure(400)

var_name = 'PTRtave01';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz = 1; % nz = 25 for PS-100
v_plot = ptr(:,:,nz);

% % for ptracer06
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


%    
    ax = subplot(13,2,[1 3 5 7 9 11]);
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(a)'])
    ylabel('Latitude')
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -38 0 38 90],...
        'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
    view(2)
    set(gca,'ColorScale','log')

% 
    ax = subplot(13,2,[17 19 21 23]);
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -1000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(b)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ylabel('Depth (m)')
    view(2)
    set(gca,'ColorScale','log')
    box on


var_name = 'PTRtave02';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz = 1; % nz = 25 for PS-100
v_plot = ptr(:,:,nz);

% % for ptracer06
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

%
    ax = subplot(13,2,[1 3 5 7 9 11]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(c)'])
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    % Create shared colorbar for first row (ax1 and ax2)
    cb1 = colorbar('Position', [0.85, 0.558, 0.02, 0.367]); % [left, bottom, width, height]
    title(cb1, 'g m^{-3}')

% 
    ax = subplot(13,2,[17 19 21 23]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -1000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(d)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    x0  =   10;
    y0  =   10;
    width   =   1300;
    height  =   800;
    set(gcf,'position',[x0,y0,width,height])   


    % Create shared colorbar for second row (ax3 and ax4)
    cb2 = colorbar('Position', [0.85, 0.174, 0.02, 0.237]);
    title(cb2, 'g m^{-3}')
    

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

    saveas(gcf, [figuredir 'figure4.png'])
    savefig(gcf, [figuredir 'figure4.fig'])


%% Figure 5

figure(500)

var_name = 'PTRtave03';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz = 25;
v_plot = ptr(:,:,nz);

%    
    ax = subplot(13,2,[1 3 5 7 9 11]);
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(a)'])
    ylabel('Latitude')
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -38 0 38 90],...
        'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
    view(2)
    set(gca,'ColorScale','log')

% 
    ax = subplot(13,2,[17 19 21 23]);
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -1000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(b)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ylabel('Depth (m)')
    view(2)
    set(gca,'ColorScale','log')
    box on


var_name = 'PTRtave04';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz = 1;
v_plot = ptr(:,:,nz);

%
    ax = subplot(13,2,[1 3 5 7 9 11]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(c)'])
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    % Create shared colorbar for first row (ax1 and ax2)
    cb1 = colorbar('Position', [0.85, 0.558, 0.02, 0.367]); % [left, bottom, width, height]
    title(cb1, 'g m^{-3}')

% 
    ax = subplot(13,2,[17 19 21 23]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -1000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(d)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    x0  =   10;
    y0  =   10;
    width   =   1300;
    height  =   800;
    set(gcf,'position',[x0,y0,width,height])   


    % Create shared colorbar for second row (ax3 and ax4)
    cb2 = colorbar('Position', [0.85, 0.174, 0.02, 0.237]);
    title(cb2, 'g m^{-3}')
    

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

    saveas(gcf, [figuredir 'figure5.png'])
    savefig(gcf, [figuredir 'figure5.fig'])



%% Figure 6

figure(600)

var_name = 'PTRtave06';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

% for ptracer06
v_plot = zeros(90,1170);
for i = 1 : 90
    for j = 1 : 1170
        kplot = 1;
        for k = 1 : 49
            if ~maskc(i,j,kplot)
                kplot = kplot;
            elseif maskc(i,j,kplot) && ~maskc(i,j,kplot+1)
                kplot = kplot;
            elseif maskc(i,j,kplot) && maskc(i,j,kplot+1)
                kplot = kplot + 1;
            end
        end
        v_plot(i,j) = ptr(i,j,kplot);
    end
end

%    
    ax = subplot(13,2,[1 3 5 7 9 11]);
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(a)'])
    ylabel('Latitude')
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-90 -38 0 38 90],...
        'yticklabel',{['90' char(176) 'S'], ['38' char(176) 'S'], ['0' char(176)], ['38' char(176) 'N'], ['90' char(176) 'N']})
    view(2)
    set(gca,'ColorScale','log')

% 
    ax = subplot(13,2,[17 19 21 23]);
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -6000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(b)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ylabel('Depth (m)')
    view(2)
    set(gca,'ColorScale','log')
    box on


var_name = 'PTRtave07';

ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

nz = 1;
v_plot = ptr(:,:,nz);

%
    ax = subplot(13,2,[1 3 5 7 9 11]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    h = trisurf(tri_z,xc_ecco,yc_ecco, v_plot);
    axis([0 360 -90 90])
    clim(color_ends_z)
    title(Tit_z)
    lighting phong
    shading interp
    colorMap = load('batlow.mat');
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(c)'])
    yline(-38,'-.w','LineWidth',1)
    box on
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    % Create shared colorbar for first row (ax1 and ax2)
    cb1 = colorbar('Position', [0.85, 0.558, 0.02, 0.367]); % [left, bottom, width, height]
    title(cb1, 'g m^{-3}')

% 
    ax = subplot(13,2,[17 19 21 23]+1);
    ax.Position(1) = ax.Position(1) - 0.07;
    hp = trisurf(tri_y,xt,zt,ptr(pplease));
    title(Tit_y)
    clim(color_ends_y);
    axis([0 360 -6000 0])
    lighting phong
    shading interp
    colormap(colorMap.batlow)
    xlabel(['Longitude' newline newline '(d)'])
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    ax.YTick = [];
    view(2)
    set(gca,'ColorScale','log')

    x0  =   10;
    y0  =   10;
    width   =   1300;
    height  =   800;
    set(gcf,'position',[x0,y0,width,height])   


    % Create shared colorbar for second row (ax3 and ax4)
    cb2 = colorbar('Position', [0.85, 0.174, 0.02, 0.237]);
    title(cb2, 'g m^{-3}')
    

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
    set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

    saveas(gcf, [figuredir 'figure6.png'])
    savefig(gcf, [figuredir 'figure6.fig'])
