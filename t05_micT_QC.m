

%8/ Apply Quality Control to MicT Measurements
%
% This script applies any quality control measures you've determined. To
% figure out which quality control parameters you want to use, play around
% with the script choose_QC_micT.m or something similar. 
%
% Input: Glider.mat, and MicT_preQC.mat from t03_micT_processing
% Output: MicT.mat with QC applied
%
% Dependencies:
%   imageQCmicT.m
%   plotQCmicT.m
%
% 2018-07-31 Updated (loosened) some of the QC controls on chi to reflect
% the new chi calculations


clear 
close all

%% %%%%%%%%%%%%%%% TOP MATTER %%%%%%%%%%%%%%%%%%%

load ../03_shear_glider/Glider.mat 

load MicT_preQC02.mat

% outfilename = 'MicT.mat'; %output filename, if saveFile = true
outfilename = 'MicT02.mat';

saveFile = false; %flag. Save results to MAT file?

nCasts = length(Glider);
epsLims = [-14 -7]; %epsilon limits when drawing figures
chiLims = [-12.5 -6]; %chi limits when drawing figures




%% %%%%%%%%%%%%%%% DEFINE QC PARAMETERS %%%%%%%%%%%%%%%%%

%GEOTRACES 2015
%
%QC01: dU/dt of glider is too big (use 90th percentile)
%QC02: within dz=15 m of top or bottom of cast
%QC03: probes disagree by more than a factor of 10
%QC04: correction term to chi is greater than 50%
%QC05: number of points in spectrum is less than n
%QC06: *REMOVED* MAD (goodness of fit of Batchelor curve) is greater than 2
%QC07: uRatio is less than 5
%
%note MAD criterion chosen heuristically, only to remove spurious
%measurements

%100 figures = epsilon figures, visualize QC
%200 figures = chi figures, visualize QC


%% Sections of epsilon and chi
%Draw the raw data as a baseline
figure(1); clf; plotQCmicT(MicTSpectra, 'eps', epsLims);
figure(2); clf; plotQCmicT(MicTSpectra, 'chi', chiLims);
drawnow

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

%% uRatio, dU/dt, and relative size of correction terms

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

%% QC01: dU/dt too big. Use a cutoff percentile

%Choose percentile above which to cutoff data
prcntCutoff=90; %make sure this is the same as for shear QC
alldat=cell2mat({QC.dUdT}'); %the non-interpolated data
prctU = prctile(alldat,prcntCutoff); %dU/dt at specified percentile
fprintf('percentile of dU/dt: %d \nvalue at percentile: %.2f m/s/day\n',...
    prcntCutoff,prctU);

%eps and chi with this QC condition applied
[QC.epsQC01] = deal(QC.eps);
[QC.chiQC01] = deal(QC.chi);
for ii=1:nCasts
    QC(ii).epsQC01(QC(ii).dUdT>prctU,:) = nan;
    QC(ii).chiQC01(QC(ii).dUdT>prctU,:) = nan;
end

%Plot the results of this QC condition
figure(101); clf; imageQCmicT(QC,'epsQC01','eps',epsLims);
figure(201); clf; imageQCmicT(QC, 'chiQC01', 'chi',chiLims);

%print
nMeas00 = sum(sum(isfinite(cell2mat({QC.eps}'))));
nMeas01 = sum(sum(isfinite(cell2mat({QC.epsQC01}'))));
fprintf('\nN epsilon at QC01: %d\n', nMeas01)
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas01, ...
    (nMeas00-nMeas01)/nMeas00*100);


drawnow
clear alldat prctU prcntCutoff

%% QC02: data too near turnaround points

%cutoff value (m)
dz=15; 

%eps and chi with this QC condition applied
[QC.epsQC02] = deal(QC.eps);
[QC.chiQC02] = deal(QC.chi);
for ii=1:nCasts
    iP = (QC(ii).P-min(QC(ii).P) < dz) | (max(QC(ii).P)-QC(ii).P < dz);
    QC(ii).epsQC02(iP,:)=nan;
    QC(ii).chiQC02(iP,:)=nan;
end

%plot
figure(102); clf; imageQCmicT(QC, 'epsQC02', 'eps', epsLims);
figure(202); clf; imageQCmicT(QC, 'chiQC02', 'chi', chiLims);

nMeas02 = sum(sum(isfinite(cell2mat({QC.epsQC02}'))));
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas02, ...
    (nMeas00-nMeas02)/nMeas00*100);

drawnow
clear dz iP

%% QC03: Two probes must be within a factor of 10, else remove both

%eps and chi with this QC condition applied
[QC.epsQC03] = deal(QC.eps);
[QC.chiQC03] = deal(QC.chi);

%if they differ by more than x10, remove both estimates
for ii=1:nCasts
    %first for chi
    chiRatio = QC(ii).chi(:,1)./QC(ii).chi(:,2);
    ibad = chiRatio>10 | chiRatio<1/10;
    QC(ii).chiQC03(ibad,:)=nan; 
    %then for epsilon
    epsRatio = QC(ii).eps(:,1)./QC(ii).eps(:,2);
    ibad = epsRatio>10 | epsRatio<1/10; 
    QC(ii).epsQC03(ibad,:)=nan;
end

%plot
figure(103); clf; imageQCmicT(QC, 'epsQC03', 'eps', epsLims);
figure(203); clf; imageQCmicT(QC, 'chiQC03', 'chi', chiLims);

%print
nMeas03 = sum(sum(isfinite(cell2mat({QC.epsQC03}'))));
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas03, ...
    (nMeas00-nMeas03)/nMeas00*100);

drawnow
clear chiRatio ibad

%% QC04: relative size of the correction terms must be <50%
%2018-07-31 - no longer apply to chi

corrCutoff=0.5;
%eps and chi with this QC condition applied
[QC.epsQC04] = deal(QC.eps);
[QC.chiQC04] = deal(QC.chi);
for ii=1:nCasts
    %apply condition
    ibad = QC(ii).corr_tot>corrCutoff;
%     QC(ii).chiQC04(ibad)=nan;
    QC(ii).epsQC04(ibad)=nan;
end

%plot
figure(104); clf; imageQCmicT(QC,'epsQC04','eps',epsLims);
figure(204); clf; imageQCmicT(QC,'chiQC04','chi',chiLims);

%print
nMeas04 = sum(sum(isfinite(cell2mat({QC.epsQC04}'))));
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas04, ...
    (nMeas00-nMeas04)/nMeas00*100);

drawnow
clear corrCutoff ibad

%% QC05: Spectrum must include at least n points
%2018-07-31 - no longer apply to chi

nCutoff=6;
%eps and chi with this QC condition applied
[QC.epsQC05] = deal(QC.eps);
[QC.chiQC05] = deal(QC.chi);
for ii=1:nCasts
    %apply condition
    ibad = MicTSpectra(ii).num_k<nCutoff;
%     QC(ii).chiQC05(ibad)=nan;
    QC(ii).epsQC05(ibad)=nan;
end

%plot
figure(105); clf; imageQCmicT(QC,'epsQC05','eps',epsLims);
figure(205); clf; imageQCmicT(QC,'chiQC05','chi',chiLims);

%print
nMeas05 = sum(sum(isfinite(cell2mat({QC.epsQC05}'))));
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas05, ...
    (nMeas00-nMeas05)/nMeas00*100);

drawnow
clear nCutoff ibad

%% QC06: MAD can't be too big

%*NOT USING THIS ANYMORE*

madCutoff = 2;
%eps and chi with this QC condition applied
[QC.epsQC06] = deal(QC.eps);
[QC.chiQC06] = deal(QC.chi);
for ii=1:nCasts
    %apply condition - only applies to epsilon
    ibad = MicTSpectra(ii).MAD>madCutoff;
    QC(ii).epsQC06(ibad)=nan;
end

%plot
figure(106); clf; imageQCmicT(QC,'epsQC06','eps',epsLims);
figure(206); clf; imageQCmicT(QC,'chiQC06','chi',chiLims);

drawnow
clear madCutoff ibad

%% QC07: uRatio is less than 5

%2018-07-31 - no longer apply to chi

uRatioCutoff = 5;
%epsilon and chi with this condition applied
[QC.epsQC07] = deal(QC.eps);
[QC.chiQC07] = deal(QC.chi);
for ii=1:nCasts
    for jj=1:2
        QC(ii).epsQC07(QC(ii).uRatio(:,jj)<=uRatioCutoff,jj) = nan;
%         QC(ii).chiQC07(QC(ii).uRatio(:,jj)<=uRatioCutoff,jj) = nan;
    end
end

%and plot
figure(107); clf; imageQCmicT(QC, 'epsQC07', 'eps', epsLims);
figure(207); clf; imageQCmicT(QC, 'chiQC07', 'chi', chiLims);

%print
nMeas07 = sum(sum(isfinite(cell2mat({QC.epsQC07}'))));
fprintf('Flagged %d (%.1f%%)\n', nMeas00-nMeas07, ...
    (nMeas00-nMeas07)/nMeas00*100);

drawnow
clear uRatioCutoff ii jj


%% Apply QC controls

%our new output structure
MicT = MicTSpectra;

%QC01: dU/dt of glider is too big
%QC02: within dz=7 m of top or bottom of cast
%QC03: probes disagree by more than a factor of 10
%QC04: correction term to chi is greater than 50% - *epsilon only*
%QC05: number of points in spectrum is less than 5 - *epsilon only*
%QC06: *REMOVED* MAD (goodness of fit of Batchelor curve) is less than 2
%QC07: uRatio is less than 5 - *epsilon only*

for ii=1:nCasts
    
    %Chi first - setup the logical conditions based on QC
    ibad = isnan(QC(ii).chiQC01) | ...  %QC01
           isnan(QC(ii).chiQC02) | ...  %QC02
           isnan(QC(ii).chiQC03) | ...  %QC03
           isnan(QC(ii).chiQC04) | ...  %QC04
           isnan(QC(ii).chiQC05) | ...  %QC05
           ... isnan(QC(ii).chiQC06) | ...  %QC06
           isnan(QC(ii).chiQC07) ;      %QC07
    %Apply the QC conditions to chi
    MicT(ii).chi(ibad) = nan;
    
    %Then Epsilon - setup the logical conditions based on QC
    ibad = isnan(QC(ii).epsQC01) | ...  %QC01
           isnan(QC(ii).epsQC02) | ...  %QC02
           isnan(QC(ii).epsQC03) | ...  %QC03
           isnan(QC(ii).epsQC04) | ...  %QC04
           isnan(QC(ii).epsQC05) | ...  %QC05
           ... isnan(QC(ii).epsQC06) | ...  %QC06
           isnan(QC(ii).epsQC07) ; ...  %QC07
    %Apply the QC conditions to epsilon
    MicT(ii).eps(ibad) = nan;
    
end

%Total number removed
nMeas = sum(sum(isfinite(cell2mat({MicTSpectra.eps}'))));
nMeasQC = sum(sum(isfinite(cell2mat({MicT.eps}'))));
fprintf('\nQC removed %d/%d (%.1f%%) of uncorrupted data\n',...
    (nMeas-nMeasQC), nMeas, (nMeas-nMeasQC)/nMeas*100)


%And draw the figures
figure(3); clf; plotQCmicT(MicT, 'eps', epsLims);
figure(4); clf; plotQCmicT(MicT, 'chi', chiLims);
drawnow

%% Save output

if saveFile
    disp('Saving...');
   save(outfilename, 'MicT'); 
end
































