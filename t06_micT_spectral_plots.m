
%The purpose of this script is to provide a detailed look at the
%dissipation estimates from a single cast, including all the spectra

%dependencies:
%   m_map
%   

%%

clear
close all
load MicT.mat
load ../03_shear_glider/Glider.mat
nCasts=length(MicT); %number of casts
jj=25; %choose a cast number. Should be in the range 1:nCasts
cm = get(groot,'DefaultAxesColorOrder');

%% make map with cast highlighted and date

%load
ibcao = load('ibcao4.mat'); %topography from IBCAO4
%setup figure
figure(31); clf; hold on; box on; grid on; zoom on
set(gca,'layer','top')
%Choose domain and the contour levels you want to draw
lolat=70.5; uplat=72.25;
lolon=-130; uplon=-123;
cont=[0 -50 -150 -250 -350 -450 -650 -1000:-500:-4000];
%extract the right part of the map
ilon = ibcao.lone>lolon & ibcao.lone<uplon;
ilat = ibcao.late>lolat & ibcao.late<uplat;
maplon = ibcao.lone(ilon);
maplat = ibcao.late(ilat);
mapZdata = ibcao.ze(ilat,ilon);
%make the map
m_proj('lambert','long',[lolon uplon],'lat',[lolat uplat]);
m_contourf(maplon,maplat,mapZdata,cont);
m_contour(maplon,maplat,mapZdata, [0 0], 'linewidth', 5, 'color', 'k')
m_grid('box','fancy','tickdir','in');
m_plot([Glider.lon],[Glider.lat],'k','linewidth',3)
m_plot(Glider(jj).lon,Glider(jj).lat,'ro','linewidth',3,'markersize',15)
m_plot(Glider(jj).lon,Glider(jj).lat,'rx','linewidth',2,'markersize',7)
%colormaps and stuff
colormap('bone'); caxis([-6000 0])
cbh=colorbar; cbh.Limits=[-4000 0];
str=datestr(datetime(datevec(mean(Glider(jj).mtime))));
m_text(Glider(jj).lon,Glider(jj).lat-0.25,str,'color','r','fontsize',15)

clear lolat lolon cont ilon ilat maplon maplat mapZdata cbh str 
clear uplat uplon ax ibcao

%% A couple of vertical profiles

figure(32); clf; ax=[];
%density
ax(1)=subtightplot(1,4,1,0.02,[0.1 0.05],[0.06 0.03]); 
hold on; box on; grid on; axis ij; zoom on; 
ylim([0 max(Glider(jj).P)+5]); 
ylabel('P (dbar)'); xlabel('\sigma (kg/m^3)');
plot(Glider(jj).rho-1000,Glider(jj).P,'linewidth',3); 
%glider speed
ax(2)=subtightplot(1,4,2,0.02,[0.1 0.05],[0.06 0.03]); 
hold on; box on; grid on; axis ij; zoom on;
ylim([0 max(Glider(jj).P)+5]); 
set(gca,'yticklabel',[]);
xlabel('U (m/s)');
plot(Glider(jj).U,Glider(jj).P,'linewidth',3); 
%epsilon 1 and 2 profiles
ax(3)=subtightplot(1,4,3,0.02,[0.1 0.05],[0.06 0.03]); 
hold on; box on; grid on; axis ij; zoom on;
ylim([0 max(Glider(jj).P)+5]); xlim(10.^[-15 -5])
set(gca,'yticklabel',[],'xscale','log');
xlabel('\epsilon (W/kg)');
plot(MicT(jj).eps(:,1), MicT(jj).P, '-s','linewidth',2,'markersize',10); 
plot(MicT(jj).eps(:,2), MicT(jj).P, '-s','linewidth',2,'markersize',10); 
%ratio of epsilon probes
ax(4)=subtightplot(1,4,4,0.02,[0.1 0.05],[0.06 0.03]); 
hold on; box on; grid on; axis ij; zoom on;
ylim([0 max(Glider(jj).P)+5]); 
set(gca,'yticklabel',[],'xscale','log');
xlabel('\epsilon_1/\epsilon_2');
plot(MicT(jj).eps(:,1)./MicT(jj).eps(:,2), MicT(jj).P, ...
    '-s','linewidth',2,'markersize',10); 

linkaxes(ax,'y');

 %% Draw Shear Spectra
 
 %Draw one figure containing all the shear spectra for this cast. To make
 %this figure legible, consecutive spectra are vertically offset from each
 %other by a constant. Both shear spectra, before and after cleaning with
 %the Goodman algorithm are shown. Also shown for each spectrum is the
 %Nasmyth spectrum for the calculated epsilon, and the 5e-11 level Nasmyth
 %spectrum for reference. Red numbers with each spectrum indicate 
 %"depth (dbar) | epsilon (W/kg)". Look at individual spectra by zooming in
 %on a few, and then using the pan tool to move up or down.
 
%Give it a minute or two...it's a big figure

%Setup figure
figure(33); clf;
npts=size(MicT(jj).eps,1); %number of spectra/epsilons
cmeps = cmocean('thermal'); %colormap
cmsc = linspace(-14,-8,length(cmeps)); %colormap scale

%setup axes for probe 1
ax2(1)=subtightplot(1,2,1,0.03,0.07,[0.07 0.05]); hold on; box on;
set(gca,'fontsize',12,'xscale','log','yscale','linear','xgrid','on');
cbh=colorbar('southoutside');
cbh.Label.String='epsilon';
colormap(cmocean('thermal'))
caxis([-14 -8])

%setup axes for probe 2
ax2(2)=subtightplot(1,2,2,0.03,0.07,[0.07 0.05]); hold on; box on;
set(gca,'fontsize',12,'xscale','log','yscale','linear','xgrid','on');
cbh=colorbar('southoutside');
cbh.Label.String='epsilon';
colormap(cmocean('thermal'))
caxis([-14 -8])

%if downcasts, swap the plotting order so we have same colour scheme throughout
if Glider(jj).dncast
    plotOrder=npts:-1:1;
else
    plotOrder=1:npts;
end

%loop through all measurements in cast
linecount=0; %counter
for hh=plotOrder
    
    linecount=linecount+1; %counter
    yshift=3.5*linecount; %vertical offset
    
    %squeeze, log, then shift the MicT and Batchelor spectra = "scaled"
    %variables. Vertical shift allows us to plot all spectra together
    Scaled.k = MicT(jj).k(hh,:);
    Scaled.Batch1 = log10(squeeze(MicT(jj).spec_batch(hh,1,:))) + yshift;
    Scaled.Batch2 = log10(squeeze(MicT(jj).spec_batch(hh,2,:))) + yshift;
    Scaled.Spec1 = log10(squeeze(MicT(jj).spec_measu(hh,1,:))) + yshift;
    Scaled.Spec2 = log10(squeeze(MicT(jj).spec_measu(hh,2,:))) + yshift;
    Scaled.Noise1 = log10(squeeze(MicT(jj).spec_noise(hh,1,:))) + yshift;
    Scaled.Noise2 = log10(squeeze(MicT(jj).spec_noise(hh,2,:))) + yshift;
    
    %Remove the long tail on the Batchelor spectra
    ibad = Scaled.Batch1<min(Scaled.Spec1);
    Scaled.Batch1(ibad)=nan;
    ibad = Scaled.Batch2<min(Scaled.Spec2);
    Scaled.Batch2(ibad)=nan;
    
    %Plot probe 2 Ref Nasmyth, Fitted Nasmyth, Raw Shear, and Cleaned Shear
    axes(ax2(1));
    icm = closeto(cmsc,log10(MicT(jj).eps(hh,1))); %colormap index
    if isnan(icm), icm=1; end;
    ph(2)=plot(Scaled.k, Scaled.Batch1, 'color', 0.5*[1 1 1], 'linewidth', 1);
    ph(3)=plot(Scaled.k, Scaled.Noise1, '--', 'color', 'k', 'linewidth', 0.5);
    ph(4)=plot(Scaled.k, Scaled.Spec1, 'color', cmeps(icm,:), 'linewidth', 2);
    %Add upper and lower integration limit
    ikmax = closeto(Scaled.k, MicT(jj).kmax(hh,1));
    ikmin = closeto(Scaled.k, MicT(jj).kmin(hh,1));
    if isfinite(ikmax)
        ph(5)=plot(MicT(jj).kmax(hh), Scaled.Spec1(ikmax), ...
            'o','color',cm(5,:), 'linewidth', 2);
        ph(6)=plot(MicT(jj).kmin(hh), Scaled.Spec1(ikmin), ...
            'o','color',cm(1,:), 'linewidth', 2);
    end
    %Add water depth and epsilon
    str=sprintf('%.1f | %.1e', MicT(jj).P(hh), MicT(jj).eps(hh,1));
    text(1.1*Scaled.k(end), Scaled.Spec1(end), str, 'color', 'r', 'fontsize', 13);
    xlim([5e-1 5e3]); pan yon; zoom yon;

    
    %Plot probe 2 Ref Nasmyth, Fitted Nasmyth, Raw Shear, and Cleaned Shear
    axes(ax2(2));
    icm = closeto(cmsc,log10(MicT(jj).eps(hh,2)));
    if isnan(icm), icm=1; end;
    plot(Scaled.k, Scaled.Batch2, 'color', 0.5*[1 1 1], 'linewidth', 1)
    plot(Scaled.k, Scaled.Noise2, '--','color', 'k', 'linewidth', 0.5)
    plot(Scaled.k, Scaled.Spec2, 'color', cmeps(icm,:), 'linewidth', 2)
    %Add upper integration limit
    ikmax = closeto(Scaled.k, MicT(jj).kmax(hh,1));
    ikmin = closeto(Scaled.k, MicT(jj).kmin(hh,1));
    if isfinite(ikmax)
        plot(MicT(jj).kmax(hh), Scaled.Spec2(ikmax), ...
            'o','color',cm(5,:), 'linewidth', 2);
        plot(MicT(jj).kmin(hh), Scaled.Spec2(ikmin), ...
            'o','color',cm(1,:), 'linewidth', 2);
    end
    %Add water depth and epsilon
    str=sprintf('%.1f | %.1e', MicT(jj).P(hh), MicT(jj).eps(hh,2));
    text(1.1*Scaled.k(end), Scaled.Spec2(end), str, 'color', 'r', 'fontsize', 13);
    xlim([5e-1 5e3]); pan yon; zoom yon;

end

%And label the axes
% axes(ax2(1));
% xlabel('k [cpm]');
% ylabel('Log_{10}(\phi) + \delta\phi   [s^{-2} cpm^{-1}]');
% tstring=sprintf('Cast %d Probe 1', jj);
% title(tstring)
% 
% legend(ph, 'Nasm 5e-11', 'Nasm Fit', 'Shear Raw', 'Shear Clean', 'Kmax')
% 
% axes(ax2(2));
% xlabel('k [cpm]');
% tstring=sprintf('Cast %d Probe 2', jj);
% title(tstring)

%Bring text to the top and link the axes
for ii=1:2
    ht=findobj(ax2(ii), 'Type', 'Text');
    uistack(ht, 'top')
end
linkaxes(ax2, 'y');
% disp('done')
%     
% clear ht tstring str ikmax yshift linecount Scaled hh ii plotOrder ph cm 
% clear npts jj
%     
%     
    
    
    
    
    
    
    
    
    
    
    