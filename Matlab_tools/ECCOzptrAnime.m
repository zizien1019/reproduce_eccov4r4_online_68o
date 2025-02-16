
clear
clc

load('XYZ_ecco.mat');

% specify where the tracer data are
datadir = '';
varname = 'PTRtave07';
Tit_part = 'Surface concentration of particles with \rho_p=1200 kg/m^3 and d=1 \mum';
% Tit_part = 'Concentration at 650 m depth of particles with \rho_p=1030 kg/m^3 and d=100 \mum';
% Tit_part = 'Surface concentration of neutrally buoyant particles';
% Tit_part = 'Surface concentration of 2D surface particles';



%% movie
% specify where to store the output movies
savemp4_dir = '';
v = VideoWriter([savemp4_dir varname],'MPEG-4');

itr = 730: 730: 227760;
% itr = 730: 730: 191260;

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;

nz=1;
tri = delaunay(xc_ecco,yc_ecco);

M(length(itr)) = struct('cdata',[],'colormap',[]);
hf = figure();

sectionL = 30;
sections = ceil(length(itr) / sectionL);


for section = 1 : sections

    if section == sections
        itr_section = itr((section-1)*sectionL+1 : length(itr));
        ptr = rdmds([datadir varname], itr_section);
    else
        itr_section = itr((section-1)*sectionL+1 : section*sectionL);
        ptr = rdmds([datadir varname], itr_section);
    end

    for nstep = 1 : length(itr_section)

    DaysAfterStart = floor(itr_section(nstep)/24);
    DateArray = datetime(727564+DaysAfterStart, 'Format','yyyy-MM-dd', 'ConvertFrom', 'datenum');
    DateString = char(DateArray);

    plotnow = ptr(:,:,:,nstep);
    plotnow(~maskc) = nan;

    h = trisurf(tri,xc_ecco,yc_ecco, plotnow(:,:,nz));
    axis([min(xc_ecco,[],'all'), max(xc_ecco,[],'all'), ...
        min(yc_ecco,[],'all'), max(yc_ecco,[],'all')])

    Tit = [Tit_part ...
        ', in ' DateString(1:7)];
    title(Tit)
    
    % adjust the max and min properly
    clim([1e-5 10]);
    % clim([5e-4 10])
    
    lighting phong
    shading interp
    colorbar EastOutside
    colormap("turbo")
    % colormap(slanCM(143))
    title(colorbar,'g/m^3','FontSize',12);
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca,'xtick',[0 90 180 270 360],...
        'xticklabel',{['0' char(176)], ['90' char(176) 'E'], ['180' char(176)], ['90' char(176) 'W'], ['0' char(176)]})
    set(gca,'ytick',[-50 0 50],...
        'yticklabel',{['50' char(176) 'S'], ['0' char(176)], ['50' char(176) 'N']})

    view(2)
    set(gca,'ColorScale','log')
    % set(gca,'ColorScale','linear')

    x0=10;
    y0=10;
    width=650;
    height=350;
    set(gcf,'position',[x0,y0,width,height])
    
    M((section-1)*sectionL+nstep) = getframe(hf);
    end

    clear ptr
end


%% save video
% v - video writer object
v.FrameRate = 20;
open(v)
writeVideo(v,M)
close(v)


