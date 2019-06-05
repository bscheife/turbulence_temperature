
%The purpose of this script is to provide a quick overview of the
%processed micT microstructure data

%dependencies:
%   Fields2Matrix()  B.Scheifele

%%

clear

% load MicT_preQC_qB2.mat
% load MicT.mat
load MicT_preQC02.mat

% MicTSpectra=MicT;

load ../03_shear_glider/Glider.mat


nCasts=length(MicTSpectra); %number of casts
cm=get(groot,'DefaultAxesColorOrder'); %colormap
set(groot,'defaultAxesTickLabelInterpreter','none'); %don't want latex
set(groot,'DefaultTextInterpreter','none'); %don't want latex



%% image as matrices

%Common variables for the axes
n = 1:nCasts;
trackDist = [Glider.trackDist];
mtime = datetime(datevec(cellfun(@mean,{MicTSpectra.mtime})));
pmax = max(cell2mat({MicTSpectra.P}'));
p = [0:0.1:pmax]; 

%The fieldnames of MicTSpectra we want to plot
fldnames1 = {'U'; 'nu'}; %one channel
fldnames2log = {'MAD'; 'eps'; 'chi'; 'X_obs'; 'X_lw'; 'X_hw'}; %2chanLogScl
fldnames2lin = {'kmin'; 'kmax'; 'kB'; 'kstar'; 'num_k'}; %two chanls
xAxis = n; %x-axis variable: select n or trackDist or mtime

figcount = 0;

%loop over 1-channel figures
for ii=1:length(fldnames1)
    figcount = figcount+1;
    figure(figcount); clf;
    %Section, top left
    ax1=subplot(4,4,[1:3 5:7 9:11]); ax1.XLabel=[];
    box on; grid on; axis ij; hold on; zoom on;
    dataMat = Fields2Matrix(MicTSpectra, fldnames1{ii}, 'P', p);
    pcolor(xAxis,p,dataMat); shading flat; colorbar; 
    ylim([0 pmax+5]); xlim([min(xAxis) max(xAxis)]);
    set(gca,'color',0*[1 1 1]); colorbar; %caxis([0.1 0.4]);
    title(fldnames1{ii}); ylabel('p dbar')
    %Mean profile, top right
    ax2=subplot(4,4,[4:4:12]); box on; grid on; hold on; zoom on; axis ij;
    plot(mean(dataMat,2,'omitnan'),p,'.','linewidth',2);
    ylim([0 pmax+5]); xlim([0.9*min(dataMat(:)) 1.1*max(dataMat(:))])
    title(fldnames1{ii});
    %Mean time series
    ax3=subplot(4,4,13:15); box on; grid on; hold on; zoom on;
    plot(xAxis, mean(dataMat,'omitnan'),'.', 'linewidth', 2)
    xlim([min(xAxis) max(xAxis)]);
    ax3.Position(3)=ax1.Position(3);
    ylabel(fldnames1{ii});
    %Histogram
    ax4=subplot(4,4,16); box on; grid on; hold on; zoom on;
    histogram(dataMat,25); 
end

%loop over 2-channel figures
for ii=1:length(fldnames2lin)
    figcount = figcount+1;
    figure(figcount); clf;
    [dataMat1 dataMat2] = Fields2Matrix(MicTSpectra, fldnames2lin{ii}, ...
        'P', p);
    dataCell = {dataMat1 dataMat2};
    for jj=1:2
        dataMat = dataCell{jj};
        %Section, top left
        ax1=subplot(8,4,(jj-1)*16+[1:3 5:7 9:11]); ax1.XTickLabel=[];
        box on; grid on; axis ij; hold on; zoom on;
        pcolor(xAxis,p,dataMat); shading flat; colorbar;
        ylim([0 pmax+5]); xlim([min(xAxis) max(xAxis)]);
        set(gca,'color',0*[1 1 1]); colorbar; %caxis([0.1 0.4]);
        if jj==1, title(fldnames2lin{ii}), end; ylabel('p dbar')
        %Mean profile, top right
        ax2=subplot(8,4,(jj-1)*16+[4:4:12]); 
        box on; grid on; hold on; zoom on; axis ij;
        plot(mean(dataMat,2,'omitnan'),p,'.','linewidth',2);
        ylim([0 pmax+5]); xlim([0.9*min(dataMat(:)) 1.1*max(dataMat(:))])
        if jj==1, title(fldnames2lin{ii}), end;
        %Mean time series
        ax3=subplot(8,4,(jj-1)*16+[13:15]); 
        box on; grid on; hold on; zoom on;
        plot(xAxis, mean(dataMat,'omitnan'),'.', 'linewidth', 2)
        xlim([min(xAxis) max(xAxis)]);
        ax3.Position(3)=ax1.Position(3);
        ylabel(fldnames2lin{ii});
        %Histogram
        ax4=subplot(8,4,(jj-1)*16+16); box on; grid on; hold on; zoom on;
        histogram(dataMat,25);
    end
end

%loop over 2-channel figures with log axes
for ii=1:length(fldnames2log)
    figcount = figcount+1;
    figure(figcount); clf;
    [dataMat1 dataMat2] = Fields2Matrix(MicTSpectra, fldnames2log{ii}, ...
        'P', p);
    dataCell = {dataMat1 dataMat2};
    for jj=1:2
        dataMat = log10(dataCell{jj});
        %Section, top left
        ax1=subplot(8,4,(jj-1)*16+[1:3 5:7 9:11]); ax1.XTickLabel=[];
        box on; grid on; axis ij; hold on; zoom on;
        pcolor(xAxis,p,dataMat); shading flat; colorbar;
        ylim([0 pmax+5]); xlim([min(xAxis) max(xAxis)]);
        set(gca,'color',0*[1 1 1]); colorbar; %caxis([0.1 0.4]);
        if jj==1, title(['log10 ' fldnames2log{ii}]), end; ylabel('p dbar')
        %Mean profile, top right
        ax2=subplot(8,4,(jj-1)*16+[4:4:12]); 
        box on; grid on; hold on; zoom on; axis ij;
        plot(mean(dataMat,2,'omitnan'),p,'.','linewidth',2);
        ylim([0 pmax+5]); xlim([0.9*min(dataMat(:)) 1.1*max(dataMat(:))])
        if jj==1, title(['log10 ' fldnames2log{ii}]), end;
        %Mean time series
        ax3=subplot(8,4,(jj-1)*16+[13:15]); 
        box on; grid on; hold on; zoom on;
        plot(xAxis, mean(dataMat,'omitnan'),'.', 'linewidth', 2)
        xlim([min(xAxis) max(xAxis)]);
        ax3.Position(3)=ax1.Position(3);
        ylabel(fldnames2log{ii});
        %Histogram
        ax4=subplot(8,4,(jj-1)*16+16); box on; grid on; hold on; zoom on;
        histogram(dataMat,25);
    end
end

clear ax1 ax2 ax3 ax4 dataCell dataMat dataMat1 dataMat2 

%% Custom figure mods

%Data limits for the specified figure numbers. Leave blank for default
figprops(3).dataLims = [0 10]; %kmin
figprops(4).dataLims = [0 100]; %kmax
figprops(5).dataLims = [0 200]; %kB
figprops(6).dataLims = [0 3]; % kstar
figprops(7).dataLims = [0 100]; %Num_k
figprops(8).dataLims = [-0.5 0.2]; %MAD
figprops(9).dataLims = [-13.5 -5]; %eps
figprops(10).dataLims = [-11 -8]; %chi
figprops(11).dataLims = [-12 -8]; %chi obs
figprops(12).dataLims = [-13 -10]; %chi lw
figprops(13).dataLims = [-13 -10]; %chi hw

%Scale the axes to the defined data limits
for jj=1:length(figprops)
    %if data limits are specified
    if ~isempty(figprops(jj).dataLims)
        lims = figprops(jj).dataLims;
        %choose figure and get axes in correct order
        figure(jj); ax = flipud(findall(gcf,'type','axes'));
        %color scale for ax 1 and 5
        for ii=[1 5], caxis(ax(ii),lims), end;
        %vert scale for ax 3 and 7
        for ii=[3 7], ylim(ax(ii),lims), end;
        %horz scale for ax 2 and 4 and 6 and 8
        for ii=[2 4 6 8], xlim(ax(ii),lims), end;
    end
end

clear ii jj lims ax

%% Degrees of Freedom check

Y = cell(nCasts,1);
for ii=1:nCasts
    
    %Get Measured and theoretical spectra
    S_meas = MicTSpectra(ii).spec_measu;
    S_theor = MicTSpectra(ii).spec_batch + MicTSpectra(ii).spec_noise;
    
    %Remove those wavenumbers not used in the Batchelor fit
    k = MicTSpectra(ii).k;
    kmin = MicTSpectra(ii).kmin;
    kmax = MicTSpectra(ii).kmax;
    eps = MicTSpectra(ii).eps;
    %loop through two probes
    for jj=1:2
        %loop through each data point in this cast
        for rr=1:size(kmax,1)
            %If no epsilon, remove S_meas
            if ~isfinite(eps(rr,jj))
                S_meas(rr,jj,:)=nan;
                continue
            end
            %Remove unused points in S_meas
            ibad = k(rr,:)<kmin(rr,jj) | k(rr,:)>kmax(rr,jj);
            S_meas(rr,jj,ibad)=nan;
        end
    end
    %Calculate ratio
    Y{ii} = S_meas(:) ./ S_theor(:);    
end
    
%Convert to numeric array; keep only finite values (not NaNs)
Y = cell2mat(Y); Y=Y(isfinite(Y));
dof = MicTSpectra(1).dof;

dof = 10;

figure(figcount+1); clf; hold on; box on; grid on; zoom on;
%Draw the histogram of d*S_obs/S_theor
histogram(dof*Y,0:2:200,'normalization','pdf')
xx = 0:0.01:200; 
plot(xx,chi2pdf(xx,dof));

    
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(groot,'defaultAxesTickLabelInterpreter','tex'); %don't want latex
set(groot,'DefaultTextInterpreter','tex'); %don't want latex














