clear
clc

load('XYZ_ecco.mat');
load('Geo_ecco.mat');

%% construct mask for each basin

Npac = ((yc_ecco>=20)&(yc_ecco<=35)&(xc_ecco>=150)&(xc_ecco<=230)) ...
    & maskc(:,:,1);

ExNpac = ((yc_ecco>=5)&(yc_ecco<=60)&(xc_ecco>=150)&(xc_ecco<=230)) ...
    & maskc(:,:,1);

Spac = ((yc_ecco>=-38)&(yc_ecco<=-18)&(xc_ecco>=180)&(xc_ecco<=260)) ...
    & maskc(:,:,1);

Natl = ((yc_ecco>=0)&(yc_ecco<=67)&(xc_ecco>=260)&(xc_ecco<=360)) ...
    & (~((yc_ecco>=0)&(yc_ecco<=15.5)&(xc_ecco>=260)&(xc_ecco<=278))) ...
    & maskc(:,:,1);

Satl = (yc_ecco<=0) & ((xc_ecco>=278)|(xc_ecco<=20.5)) ...
    & maskc(:,:,1);

Indi = ((xc_ecco>=20.5)&(yc_ecco<=26)&(xc_ecco<=140.5)) ...
    & (~((yc_ecco>=-6.5)&(yc_ecco<=26)&(xc_ecco>=103.5)&(xc_ecco<=140.5))) ...
    & maskc(:,:,1);

Medi = (yc_ecco>=30)&(yc_ecco<=50)&(xc_ecco>=0)&(xc_ecco<=50) ...
    & maskc(:,:,1);

Glob = maskc(:,:,1);

% look at the values in North Pacific in the following codes
here = Npac;


%%

% specify where tracer data are stored
datadir = '';
ptr_num = 1;
[ptr_3d,itr] = rdmds([datadir 'PTRtave0' num2str(ptr_num)], nan);

% computes the avaraged value over a specific basin
% specify which basin to look at in the previous section 
ptr_have = zeros(length(itr));
Asum = tensorprod(double(here), RAC, [1 2],[1 2]);

for i = 1 : length(itr)
    ptr_temp = ptr_3d(:,:,:,i);
    ptr_2d = ptr_temp(:,:,1);
    ptr_have(i) = tensorprod(ptr_2d, RAC, [1 2],[1 2]) ./ Asum;
end

xc_ecco(xc_ecco<0) = xc_ecco(xc_ecco<0) + 360;



%%
    DaysAfterStart = floor(itr/24);
    Month = mod(ceil(DaysAfterStart/(365/12)-1), 12)+1;
    Year = 1992 + DaysAfterStart / 365;
    
    yr13to17 = 253:312;
    plot(Year(yr13to17),ptr_have(yr13to17));
    set(gca,'xtick',2013:0.5:2018,...
        'xticklabel',{'','2013','','2014','','2015','','2016','','2017',''})
    xlabel('Year')
    ylabel('\tau (g/m^3)')
    % title('Monthly average concentration in North Pacific')

    x0=10;
    y0=10;
    width=480;
    height=200;
    set(gcf,'position',[x0,y0,width,height])
    % saveas(gcf, [figuredir 'check_' var_name '_z.png'])




