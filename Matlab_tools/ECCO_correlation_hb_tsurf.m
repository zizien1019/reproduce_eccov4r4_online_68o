clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

%% 1. load tracer data

% specify where the data are stored
datadir = ''; 

% monthly average tracer values
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


%% 2. load density to determine hb

% directory containing the file
% specify where the data are stored
rho_dir = '';

% Load Files
% monthly average density anomoly
% 2012 Jan - 2017 Dec
% 72 records in total
itr = [176052 176748 177492 178212 178956 179676 180420 181164 181884 182628 183348 184092 184836 185508 186252 186972 187716 188436 189180 189924 190644 191388 192108 192852 193596 194268 195012 195732 196476 197196 197940 198684 199404 200148 200868 201612 202356 203028 203772 204492 205236 205956 206700 207444 208164 208908 209628 210372 211116 211812 212556 213276 214020 214740 215484 216228 216948 217692 218412 219156 219900 220572 221316 222036 222780 223500 224244 224988 225708 226452 227172 227903];

rho = rdmds([rho_dir 'RHOAnoma_mon_mean'], itr);


hb = zeros(90, 1170, 72);

for t = 1 : 72
    for i = 1 : 90
        for j = 1 : 1170

            r = squeeze(rho(i,j,:,t));
            drdz = gradient(r) ./ gradient(zc_ecco);
            [~,index] = max(-drdz);
            hb(i,j,t) = zc_ecco(index);
        end
    end
end

clear rho

%% plot to check
% tri = delaunay(xc_ecco,yc_ecco);
% ntime = 65;
% trisurf(tri,xc_ecco,yc_ecco, MXLD(:,:,ntime))
% lighting phong
% shading interp
% colorbar EastOutside
% % clim([-500 0])
% % colormap(slanCM(104))
% view(2)

%%
% 3. compute correlation
% along the third dimension (time)!
% syntax : correlation = corr([1; 2; 3], [30; 90; 80]);
% input matrices must be in columns

cor_tzave_hb = zeros(90, 1170);
for i = 1 : 90
    for j = 1 : 1170
        tzave = squeeze(ptr_surf(i,j,:));
        hb_here = squeeze(hb(i,j,:));
        cor_tzave_hb(i,j) = corr(tzave(:), hb_here);
    end
end


%%

tri = delaunay(xc_ecco,yc_ecco);
% specify where to store the figures
figuredir = '';
figure(1)
    subplot(14,1,[2:13])
trisurf(tri,xc_ecco,yc_ecco, cor_tzave_hb)
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
    saveas(gcf, [figuredir '10mum_corr_hb_tsurf.png'])