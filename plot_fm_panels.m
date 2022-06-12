function plot_fm_panels(x,y,z,title_str,Nx,Ny,varargin);
% plot out the InSAR on panels using imagesc
% 
% plot_intf_panels_ll(x,y,z,title_str,Nx,Ny);
% availabel options:
%
%   'ColorRange'   ------  [cmin cmax] for color range
%   'Region'       ------  [xmin xmax ymin ymax] for boundary
%   'FigureName'   ------  ['name_of_figure']
%   'color_map'    ------  ['string of the GMT colormap']
%   'FontSize'     ------  [font_size]
%   'TextSize'     ------  [font size for the text]
%   'FontName'     ------  ['font_name']
%   'ColorbarTitle'------  ['title on the colorbar']
%   'Units'        ------  ['units of the data']
%   'PlotPoint'    ------  [lon_pt,lat_pt]  can be use to plot fault trace
%   'dir_figure'   ------  ['dir_figure_out']
%


% some default values....
%
%%%%%  size of canvas %%%%%
xwidth=0.50; % units are normalized to the screen size
ylength=1.0; % units are normalized to the screen size
ystretch=4;

xx=x{1};
yy=y{1};
xmin=min(xx);
xmax=max(xx);
ymin=min(yy);
ymax=max(yy);

region1=[xmin xmax ymin ymax];

color_flag=0;

figure_name=['focmec_panel'];
fontsize=8;
fontname=['Helvetica'];
tstr=['fault segments']; % title str
units=['mm'];
color_map=['jet'];
plot_pt=[NaN,NaN];
fontsize_text=10;
dir_figure=['.'];
region2=[];
if (~isempty(varargin))
  for c=1:floor(length(varargin)/2);
     
        switch varargin{c*2-1};
            case 'ColorRange'
                color_range=real(varargin{c*2});
                cmin=color_range(1);
                cmax=color_range(2);
                color_flag=1;
            case 'Colormap'
                color_map=char(varargin{c*2});
            case 'Region'
                region2=real(varargin{c*2});
            case 'FigureName'
                figure_name=char(varargin{c*2});
            case 'FontSize'
                fontsize=real(varargin{c*2});
            case 'TextSize'
                fontsize_text=real(varargin{c*2});
            case 'FontName'
                fontname=char(varargin{c*2});
            case 'ColorbarTitle'
                tstr=char(varargin{c*2});
            case 'Units'
                units=char(varargin{c*2});
            case 'PlotPoint'
                plot_pt=real(varargin{c*2});
            case 'dir_figure'
                dir_figure=char(varargin{c*2});
            case 'aspect_ratio'
                aspect_ratio=char(varargin{c*2});
            otherwise
                 error(['Unrecognized Keywords: ',varargin{c*2-1}]);
        end

  end
end

if length(region2)<4;
  region = region1;
else
  region = region2(1:4);
end

xmin=region(1);
xmax=region(2);
ymin=region(3);
ymax=region(4);

if strcmp(dir_figure(length(dir_figure)),'/');
   dir_figure=dir_figure(1:length(dir_figure)-1);
end

if ~exist(dir_figure,'dir');
  mkdir (dir_figure);
end

title_colorbar=[tstr,' (',units,')'];

set(0,'defaultAxesFontSize',fontsize)
set(0,'defaultAxesFontName',fontname)

if (color_flag<0);
 % determine the default color range
 ndata=length(z);
 zmin=zeros(ndata,1);
 zmax=zeros(ndata,1);
 for k=1:ndata;
   znow=z{k};
   zgood=znow(~isnan(znow));
   zmin_now=min(zgood);
   zmax_now=max(zgood);
   zmin(k)=zmin_now;
   zmax(k)=zmax_now;
end
  cmin_tmp=min(zmin);
  cmax_tmp=max(zmax);
  rz=cmax_tmp-cmin_tmp;
  cmin=cmin_tmp+0.05*rz;
  cmax=cmax_tmp-0.05*rz;
end

xrng=region(2)-region(1);
yrng=region(4)-region(3);

ratio_xy=yrng/xrng;  %ratio of y-axis to x-axis
yunit=0.8/Ny;
xunit=yunit/ratio_xy;


ysize=yunit*Ny;
xsize=xunit*Nx;
if (xsize > 0.8);
   xsize= 0.8;
   xunit=xsize/Nx
   yunit=xunit*ratio_xy
end

yunit=yunit*ystretch;

%data_flt=load('fault_new.txt');
lon_flt=plot_pt(:,1);
lat_flt=plot_pt(:,2);

%strcmd=['makecpt -C',color_map,' -T',num2str(cmin),'/',num2str(cmax),'/1',' -D >los.cpt'];
%disp(strcmd);
%system([strcmd]);
%cmap=importcpt('los.cpt');
cmap='jet';

Nintf=length(z);
disp(['Number of total images to display: ',num2str(Nintf)]);

disp(['Saving Figures to ',figure_name]);
Npage=Nx*Ny;
nsubplot=0;

xstart=0.5-(Nx/2)*xunit;
ystart=0.89-yunit; 
for i=1:Ny;
 for j=1:Nx;
   nsubplot=nsubplot+1;
   x1=xstart + (j-1)*xunit;
   y1=ystart - (i-1)*yunit;
   psv{nsubplot}=[x1,y1,xunit*0.98,yunit*0.98];
 end 
end

k=0;
nplot=0;

nfigure=0;
for i=1:Nintf;
  if (mod(i,Npage)==0);
    k=k+1;
    hf= figure('units','normalized','outerposition',[0.1 0 xwidth ylength]);
    set(hf,'Visible','off')
    nfigure=nfigure+1;
    for nn=1:Npage;
       mm=(k-1)*Npage+nn;
    if (mm <=Nintf)
    this_x=x{mm};
    this_y=y{mm};
    this_z=z{mm};
    xmax=max(this_x);
    xmin=min(this_x);

    ymax=max(this_y);
    ymin=min(this_y);

    Rx=xmax-xmin;
    Ry=ymax-ymin;
     subplot('position',psv{nn});
     hi=imagesc(this_x,this_y,this_z);
     
     colormap(cmap), hold on;
     caxis([cmin cmax]), hold on;
     plot(lon_flt,lat_flt,'k.'), hold on;

    set(gca,'XAxisLocation','top'), hold on;
    set(hi,'alphadata',~isnan(this_z)), hold on;

    text(xmin+0.5*Rx,ymin+0.85*Ry,title_str{mm},'interpreter','none','HorizontalAlignment','center','FontWeight','Bold','FontSize',fontsize_text), hold on;
%    axis equal, hold on;
%    axis fill, hold on
%    axis manual, hold on
    daspect([ystretch 1  1]), hold on
    axis(region);


    if (nn==1);
%      hcolor=colorbar('north','Position',[xstart+(Nx-1)*xunit 0.9,xunit,0.01]);
%      xlabel(hcolor,title_colorbar);
    end
    if (nn~=1);
      set(gca,'YTick',[])
      set(gca,'XTick',[])
    end
     set(gca,'YDir','normal')
%      text(xmin+0.5*Rx,ymin+0.9*Ry,title_str{mm},'interpreter','none','HorizontalAlignment','center','FontWeight','Bold','FontSize',20)
    end
    end
%   saveas(hf,[figure_name,'_',sprintf('%03d',nfigure)],'png');
   saveas(hf,[dir_figure,'/',figure_name,'_',sprintf('%03d',nfigure)],'png');

  end
end

mmax=k*Npage;
if (mmax<Nintf);
   hf= figure('units','normalized','outerposition',[0.1 0 xwidth ylength]);
   set(hf,'Visible','off')
    nfigure=nfigure+1;
    for kk=1:(Nintf-mmax);
    this_x=x{mmax+kk};
    this_y=y{mmax+kk};
    this_z=z{mmax+kk};
     subplot('position',psv{kk});
%    hi=imagesc(this_x,this_y,this_z);
    hold on
    colormap(cmap), hold on;
    caxis([cmin cmax]), hold on;
    plot(lon_flt,lat_flt,'k.'), hold on;
    set(gca,'XAxisLocation','top'), hold on;
    set(hi,'alphadata',~isnan(this_z)), hold on;
%    axis equal
%    axis manual
    axis(region), hold on;
   daspect([ystretch 1  1]), hold on;


    if (kk==1);
     hcolor=colorbar('north','Position',[xstart+(Nx-1)*xunit 0.9,xunit,0.01]);

    xlabel(hcolor,title_colorbar);
    end

    if (kk~=1);
     set(gca,'YTick',[])
     set(gca,'XTick',[])
    end
    set(gca,'YDir','normal')
%    xmax=max(this_x);
%    xmin=min(this_x);
%    ymax=max(this_y);
%    ymin=min(this_y);

    Rx=xmax-xmin;
    Ry=ymax-ymin;
    text(xmin+0.5*Rx,ymin+0.85*Ry,title_str{mmax+kk},'interpreter','none','HorizontalAlignment','center','FontWeight','Bold','FontSize',fontsize_text)

 end
%  saveas(hf,[figure_name,'_',sprintf('%03d',nfigure)],'png');
  saveas(hf,[dir_figure,'/',figure_name,'_',sprintf('%03d',nfigure)],'png');
end

