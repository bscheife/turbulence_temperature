


clear 
close all

%% %%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%

load ../03_shear_glider/Glider.mat 
load MicT_preQC2.mat

nCasts = length(Glider);
epsLims = [-14 -6]; %epsilon limits when drawing figures
chiLims = [-12.5 -6]; %chi limits when drawing figures

%000 figures = develop QC parameters
%100 figures = epsilon figures, visualize QC
%200 figures = chi figures, visualize QC


%% Sections of epsilon and chi
%Draw the raw data as a baseline
figure(100); clf; plotQCmicT(MicTSpectra, 'eps', epsLims);
figure(200); clf; plotQCmicT(MicTSpectra, 'chi', chiLims);

%% New structure

% Make new disposable structure that we can use to assess the Quality
% Control information we need. Should have all variables at the same
% resolution as epsilon.

QC(nCasts).P=[];
%Get measurements, 
[QC.P] = deal(MicTSpectra.P);
[QC.mtime] = deal(MicTSpectra.mtime);
[QC.U] = deal(MicTSpectra.U);
[QC.eps] = deal(MicTSpectra.eps); 
[QC.chi] = deal(MicTSpectra.chi);

%% Calculate uRatio, dU/dt, and relative size of correction terms

%loop through casts, interpolate onto epsilon grid as needed
for ii=1:nCasts
    %N2
    QC(ii).N2 = interp1(Glider(ii).Pmid,Glider(ii).N2,MicTSpectra(ii).P);
    %this isn't going to work if N2<0; keep this step in mind for later...
    QC(ii).N2(QC(ii).N2<0)=nan;
    %Turbulent velocity scale is roughly (epsilon/N)^(1/2)
    QC(ii).uTurb = sqrt(QC(ii).eps./[QC(ii).N2 QC(ii).N2]);
    %Ratio of Glider to Turbulent velocity
    QC(ii).uRatio = [QC(ii).U QC(ii).U] ./ QC(ii).uTurb;
    %Rate of change of glider speed. Have to interpolate back to right grid
    QC(ii).dUdT = abs(diff(QC(ii).U))./diff(QC(ii).mtime); %m/s/day
    QC(ii).dUdT = interp1(midpoints(QC(ii).P), QC(ii).dUdT, QC(ii).P);
    %Relative size of correction terms
    QC(ii).corr_lw = MicTSpectra(ii).X_lw ./ MicTSpectra(ii).chi;
    QC(ii).corr_hw = MicTSpectra(ii).X_hw ./ MicTSpectra(ii).chi;
    QC(ii).corr_tot = QC(ii).corr_lw + QC(ii).corr_hw;
end

%% Plot dU/dt. Use a cutoff percentile

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];
%Choose percentile above which to cutoff data
prcntCutoff=90; %make sure this is the same as for shear QC

%image dU/dt
figure(1); clf; 
dat1 = Fields2Matrix(QC,'dUdT','P',p);
subtightplot(2,2,1); hold on; box on; axis ij; zoom on; 
set(gca,'color','k'); ylim([0 max(p)])
pcolor(n,p,dat1); shading flat; caxis([0 50]);
title('dU/dt')

%Non-interpolated data
alldat=cell2mat({QC.dUdT}'); %the non-interpolated data
prctU = prctile(alldat,prcntCutoff); %dU/dt at specified percentile

%image again, showing what will be removed
subtightplot(2,2,2); hold on; box on; axis ij; zoom on; 
set(gca,'color','k','yticklabel',[]); ylim([0 max(p)])
pcolor(n,p,dat1); shading flat; caxis([prctU prctU+1]);
title('dU/dt to remove')

%draw histogram
subplot(2,1,2); zoom on; grid on
histogram(alldat,1:50)
xlabel('dU/dt (m/s/day)')

fprintf('percentile used: %d \nvalue at percentile: %.2f m/s/day\n',...
    prcntCutoff,prctU);

%image epsilon and chi with these removed
%epsQC01 : epsilon with values removed if dU/dt is too big

[QC.epsQC01] = deal(QC.eps);
[QC.chiQC01] = deal(QC.chi);
for ii=1:nCasts
    QC(ii).epsQC01(QC(ii).dUdT>prctU,:) = nan;
    QC(ii).chiQC01(QC(ii).dUdT>prctU,:) = nan;
end

figure(101); clf; imageQCmicT(QC,'epsQC01','eps',epsLims);
figure(201); clf; imageQCmicT(QC, 'chiQC01', 'chi',chiLims);

clear sh1Mat sh2Mat ax n p dat1 dat2 alldat prctU prcntCutoff


%% Plot uRatio
%Fer et al 2014 use uRatio < 20 as a cutoff. For us, uRatio<10 looks like
%it effectively removes spurious epsilons

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];

figure(2); clf; 
[dat1 dat2] = Fields2Matrix(QC,'uRatio','P',p);

subtightplot(3,2,1); hold on; box on; axis ij; zoom on;
set(gca,'color','k'); title('U/u_t');
pcolor(n,p,dat1); shading flat; colorbar; caxis([0 500]);

subtightplot(3,2,2); hold on; box on; axis ij; zoom on;
pcolor(n,p,dat1); shading flat; caxis([9.99 10]);
set(gca,'color','k', 'yticklabel',[]); colormap(flipud(parula));
title('U/u_t to remove')

subtightplot(3,2,3); hold on; box on; axis ij; zoom on;
pcolor(n,p,dat2); shading flat; colorbar; caxis([0 500]);
set(gca,'color','k'); 

subtightplot(3,2,4); hold on; box on; axis ij; zoom on;
pcolor(n,p,dat2); shading flat; caxis([9.99 10]);
set(gca,'color','k','yticklabel',[])

subplot(3,2,[5:6]); zoom on; grid on
histogram(cell2mat({QC.uRatio}'),1:2:200)
xlabel('R = U/u_t')


%calculate epsilon and chi with this condition applied
uRatioCutoff = 10;
[QC.epsQC02] = deal(QC.eps);
[QC.chiQC02] = deal(QC.chi);
for ii=1:nCasts
    for jj=1:2
        QC(ii).epsQC02(QC(ii).uRatio(:,jj)<=uRatioCutoff,jj) = nan;
        QC(ii).chiQC02(QC(ii).uRatio(:,jj)<=uRatioCutoff,jj) = nan;
    end
end

%and plot
figure(102); clf; imageQCmicT(QC, 'epsQC02', 'eps', epsLims);
figure(202); clf; imageQCmicT(QC, 'chiQC02', 'chi', chiLims);

clear n p dat1 dat2 sh1Mat sh2Mat

%% data near turnaround points
%epsQC03 : remove data closer to turnaround than dz

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];

dz=7; %if within dz meters of turnaround point, remove data
[QC.epsQC03] = deal(QC.eps);
[QC.chiQC03] = deal(QC.chi);
for ii=1:nCasts
    iP = (QC(ii).P-min(QC(ii).P) < dz) | (max(QC(ii).P)-QC(ii).P < dz);
    QC(ii).epsQC03(iP,:)=nan;
    QC(ii).chiQC03(iP,:)=nan;
end

%plot
figure(103); clf; imageQCmicT(QC, 'epsQC03', 'eps', epsLims);
figure(203); clf; imageQCmicT(QC, 'chiQC03', 'chi', chiLims);

clear n p dat1 dat2 sh1Mat sh2Mat iP dz ii jj

%% Two probes must be within a factor of 10
%epsQC04: If the two probes differ by more than a factor of 10, we keep
%only the lower of the estimates. Do this separately for epsilon and chi

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];

[QC.epsQC04] = deal(QC.eps);
[QC.chiQC04] = deal(QC.chi);

for ii=1:nCasts
    %first for chi
    chiRatio = QC(ii).chi(:,1)./QC(ii).chi(:,2);
    iR1 = chiRatio>10; %probe 1 is bigger
    iR2 = chiRatio<1/10; %probe 2 is bigger
    QC(ii).chiQC04(iR1,1)=nan;
    QC(ii).chiQC04(iR2,2)=nan;    
    %then for epsilon
    epsRatio = QC(ii).eps(:,1)./QC(ii).eps(:,2);
    iR1 = epsRatio>10; %probe 1 is bigger
    iR2 = epsRatio<1/10; %probe 2 is bigger
    QC(ii).epsQC04(iR1,1)=nan;
    QC(ii).epsQC04(iR2,2)=nan;
end

%plot
figure(104); clf; imageQCmicT(QC, 'epsQC04', 'eps', epsLims);
figure(204); clf; imageQCmicT(QC, 'chiQC04', 'chi', chiLims);

clear n p sh1Mat sh2Mat ii epsRatio iR1 iR2 chiRatio

%% Check relative size of the correction terms

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];

figure(3); clf; 
%low and high wavenumber corrections, separately
subtightplot(2,2,1,[0.07 0.01]); hold on; box on; grid on; 
histogram(cell2mat({QC.corr_lw}'),0:0.02:1);
histogram(cell2mat({QC.corr_hw}'),0:0.02:1);
legend('X_{lw}/X', 'X_{hw}/X');
set(gca,'yticklabel',[]);
%total correction, relative size
subtightplot(2,2,2,[0.07 0.01]); hold on; box on; grid on;
histogram(cell2mat({QC.corr_tot}'),0:0.02:1);
legend('X_{corr}/X')
set(gca,'yticklabel',[]);
%image correction term, both probes
[dat1 dat2] = Fields2Matrix(QC,'corr_tot','P',p);
subtightplot(2,2,3,[0.07 0.01]); hold on; box on; grid on; axis ij; zoom on;
pcolor(n,p,dat1); shading flat; set(gca,'color','k')
colorbar; 
caxis([0.49 0.50])
subtightplot(2,2,4,[0.07 0.01]); hold on; box on; grid on; axis ij; zoom on;
pcolor(n,p,dat2); shading flat; set(gca,'color','k')
colorbar; 
caxis([0.49 0.50])

%use condition if Xcorr > 50%, discard data
[QC.epsQC05] = deal(QC.eps);
[QC.chiQC05] = deal(QC.chi);
for ii=1:nCasts
    %apply condition
    ibad = QC(ii).corr_tot>0.5;
    QC(ii).chiQC05(ibad)=nan;
    QC(ii).epsQC05(ibad)=nan;
end

%plot
figure(105); clf; imageQCmicT(QC,'epsQC05','eps',epsLims);
figure(205); clf; imageQCmicT(QC,'chiQC05','chi',chiLims);

clear n p ibad

%% Number of points being integrated

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];

figure(4); clf; box on; grid on; hold on; zoom on
histogram(cell2mat({MicTSpectra.num_k}'),1:500)
xlabel('number of points in spectrum')

%Okay, looks like we want a condition for n<5
[QC.epsQC06] = deal(QC.eps);
[QC.chiQC06] = deal(QC.chi);
for ii=1:nCasts
    %apply condition
    ibad = MicTSpectra(ii).num_k<4;
    QC(ii).chiQC06(ibad)=nan;
    QC(ii).epsQC06(ibad)=nan;
end

%plot
figure(106); clf; imageQCmicT(QC,'epsQC06','eps',epsLims);
figure(206); clf; imageQCmicT(QC,'chiQC06','chi',chiLims);

clear n p ibad

%% Check MAD

n=1:nCasts;
p=[0:0.1:max(cell2mat({QC.P}'))];
[dat1 dat2] = Fields2Matrix(MicTSpectra,'MAD','P',p);

figure(5); clf; 
%histogram
subplot(211); box on; grid on; hold on; zoom on;
histogram(cell2mat({MicTSpectra.MAD}'), 10.^[-1.5:0.05:1])
xlabel('MAD'); set(gca,'xscale','log');
%sections
subtightplot(2,2,3); box on; grid on; hold on; axis ij; zoom on;
pcolor(n,p,log10(dat1)); shading flat; set(gca,'color','k')
caxis([0.99*log10(2) log10(2)])
caxis(log10([0.1 3]))
subtightplot(2,2,4); box on; grid on; hold on; axis ij; zoom on;
pcolor(n,p,log10(dat2)); shading flat; set(gca,'color','k')
caxis([0.99*log10(2) log10(2)])
caxis(log10([0.1 3]))

%Hmmm let's use a cutoff of MAD=2 for now
madCutoff = 2;
[QC.epsQC07] = deal(QC.eps);
[QC.chiQC07] = deal(QC.chi);
for ii=1:nCasts
    %apply condition
    ibad = MicTSpectra(ii).MAD>2;
    QC(ii).chiQC07(ibad)=nan;
    QC(ii).epsQC07(ibad)=nan;
end

%plot
figure(107); clf; imageQCmicT(QC,'epsQC07','eps',epsLims);
figure(207); clf; imageQCmicT(QC,'chiQC07','chi',chiLims);


clear n p dat1 dat2 madCutoff ibad


%% Apply QC controls

MicT = MicTSpectra;

%QC01 - dU/dt too big
%QC02 - U/u_t too small
%QC03 - too close to turnaround point
%QC04 - estimates differ by >x10
%QC05 - correction to chi is >50%
%QC06 - fewer than n points in spectrum
%QC07 - MAD is too big

for ii=1:nCasts
    %Chi first - setup the logical conditions based on QC
    ibad = isnan(QC(ii).chiQC01) | ...  %QC01
           isnan(QC(ii).chiQC02) | ...  %QC02
           isnan(QC(ii).chiQC03) | ...  %QC03
           isnan(QC(ii).chiQC04) | ...  %QC04
           isnan(QC(ii).chiQC05) | ...  %QC05
           isnan(QC(ii).chiQC06) | ...  %QC06
           isnan(QC(ii).chiQC07) ;      %QC07
    %Apply the QC conditions to chi
    MicT(ii).chi(ibad) = nan;
    %Then Epsilon - setup the logical conditions based on QC
    ibad = isnan(QC(ii).epsQC01) | ...  %QC01
           isnan(QC(ii).epsQC02) | ...  %QC02
           isnan(QC(ii).epsQC03) | ...  %QC03
           isnan(QC(ii).epsQC04) | ...  %QC04
           isnan(QC(ii).epsQC05) | ...  %QC05
           isnan(QC(ii).epsQC06) | ...  %QC06
           isnan(QC(ii).epsQC07) ; ...  %QC07
    %Apply the QC conditions to epsilon
    MicT(ii).eps(ibad) = nan;
    
end

%And draw the figures
figure(110); clf; plotQCmicT(MicT, 'eps', epsLims);
figure(210); clf; plotQCmicT(MicT, 'chi', chiLims);


%% Can we get away with a different uRatio cutoff

%Hmmm we're losing some of the most interesting data because of the uRatio
%cutoff. This criterion is based on Fer 2014 Fig 8, so let's also make this
%plot and see if we can get away with a smaller uRatio cutoff

MicT_temp(nCasts).eps=[];
%Get measurements, 
[MicT_temp.uRatio] = deal(QC.uRatio);
[MicT_temp.eps] = deal(MicTSpectra.eps); 

%Apply all QC conditions *EXCEPT* QC02 (and QC04)
for ii=1:nCasts
    %Epsilon - setup the logical conditions based on QC
    ibad = isnan(QC(ii).epsQC01) | ...  %QC01
           ... %isnan(QC(ii).epsQC02) | ...  %QC02
           isnan(QC(ii).epsQC03) | ...  %QC03
           isnan(QC(ii).epsQC04) | ...  %QC04
           isnan(QC(ii).epsQC05) | ...  %QC05
           isnan(QC(ii).epsQC06) | ...  %QC06
           isnan(QC(ii).epsQC07) ; ...  %QC07
    %Apply the QC conditions to epsilon
    MicT_temp(ii).eps(ibad) = nan;
end

uRatio = mean(cell2mat({MicT_temp.uRatio}'),2);
eps = cell2mat({MicT_temp.eps}');

%reorder the elements so the bigger epsilon is always in the left column
for ii=1:length(eps)
    if eps(ii,2)>eps(ii,1)
        eps(ii,:) = fliplr(eps(ii,:));
    end
end

epsRatio = eps(:,1)./eps(:,2);

[epsRatioBin, uRatioBin]=binProfile(uRatio, epsRatio, 1);

%%

figure(6); clf; box on; hold on; grid on; zoom on;
plot(uRatio, epsRatio, '.')
xlim([0 40]); ylim([0 15])
plot(uRatioBin, epsRatioBin, 'ks', 'markersize', 10, 'linewidth', 3)

title('micT'); xlabel('U/u_t'); ylabel('eps1/eps2');
%Yah okay this doesn't make a difference at all...we defn don't need that
%condition



































