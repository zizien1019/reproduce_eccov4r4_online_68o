
clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

Nstep = [219730 220460 221190 221920 222650 223380 224110 224840 225570 226300 227030 227760];

datadir = '/Users/zengzien/Code_Research/GreatLakes/small_26yrs/ptr_3d/';
var_name = 'PTRtave05';
ptr_all = rdmds([datadir var_name], Nstep);
ptr = mean(ptr_all,4,"omitmissing");
ptr(~maskc) = nan;

ptr_mean = mean(ptr(:,:,1),"all","omitmissing");

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0);

%%
tri = delaunay(xc_ecco,yc_ecco);
hf = figure();

    h = trisurf(tri,xc_ecco,yc_ecco,ptr(:,:,1)/ptr_mean);

    axis([min(xc_ecco,[],'all'), max(xc_ecco,[],'all'), ...
        min(yc_ecco,[],'all'), max(yc_ecco,[],'all')])

    % Tit = ['surface \tau (ppm), ' ...
    %     'neutrally buoyant' ...
    %     ',26 yrs'];
    % title('(a)')
    
    
    lighting phong
    shading interp
    % colorbar EastOutside
    % colormap("turbo")
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',[-180 -90 0 90 180],...
        'xticklabel',{['180' char(176)], ['90' char(176) 'W'], ['0' char(176)], ['90' char(176) 'E'], ['180' char(176)],})
    set(gca,'ytick',[-50 0 50],...
        'yticklabel',{['50' char(176) 'S'], ['0' char(176)], ['50' char(176) 'N']})

    view(2)
    set(gca,'ColorScale','log')
    % set(gca,'ColorScale','linear')

    x0=10;
    y0=10;
    width=350;
    height=260;
    set(gcf,'position',[x0,y0,width,height])
    

    colors = {'#480000', '#5e0000', '#760000', '#910000', '#ad0000', ...
        '#cb0000', '#e90000', '#ff0000', '#ff0000', '#ff0000', '#ff0012', ...
        '#ff303f', '#ff4f58', '#ff6a70', '#ff8589', '#ffa2a5', '#ffbfc1', ...
        '#c8c8c8', '#c8c8c8', '#c8c8c8', '#c8c8c8', ...
        '#c1c6ff', '#a2bdff', '#83b2ff', '#65aaff', '#47a0ff', ...
        '#2496ff', '#0084ff', '#0075ff', '#0065ff', '#0057fd', '#004ddd', ...
        '#0042c0', '#0038a4', '#002f89', '#00266f', '#001e59', '#001743'};

J = customcolormap(linspace(0,1,38), colors);
colorbar; colormap(J); 
clim([5e-3 2e2]);

%% movie

% savemp4_dir = '20240715_sizes_fine/';
% 
% % v = VideoWriter([savemp4_dir, 'monthly_' num2str(d_p(ptr_num)*1e6) 'mum_diff_colorbar'],'MPEG-4');
% v = VideoWriter([savemp4_dir, 'monthly_neu_diff_colorbar'],'MPEG-4');
% 
% % v - video writer object
% open(v)
% writeVideo(v,M)
% close(v)

%%
% 
% hp = ECCOzplotGlobal(ptr(:,:,:,end),tri,xc_ecco,yc_ecco,1,1);
%     axis([min(xc_ecco,[],'all'), max(xc_ecco,[],'all'), ...
%         min(yc_ecco,[],'all'), max(yc_ecco,[],'all')])
%     caxis([1e-6 1e-1])
%     lighting phong
%     shading interp
%     colorbar EastOutside
%     colormap("turbo")
%     xlabel('Longitude(^o)')
%     ylabel('Latitude(^o)')
%     view(2)
%     set(gca,'ColorScale','log')


