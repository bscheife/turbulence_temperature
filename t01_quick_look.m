
%A quick script to plot all the micT temperature profiles at once (with the
%glider temperature profiles). We scale down the resolution before plotting
%so things don't get too slow. You can modify the factor by which the
%plotting resolution is scaled down from the data resolution. See Section
%-Front Matter-. Profiles are offset from each other by a user-defined
%factor. Zoom in on the picture to view individual profiles.

clear 
% close all

%% Front Matter

%Input data file(s). Specify filenames as: {'fname1.mat'; 'fname2.mat'}
fnames = {'GlidMic01.mat'; 'GlidMic02.mat'; 'GlidMic03.mat'}; %file name(s)
% fnames = {'GlidMicTest.mat'};
pth = '../02_collect/'; %path to prepend to input file(s)

%Some plotting options
delT=1.5; %temperature offset between profiles
scDn=25; %scale down MR resolution by this factor, else too much data
minP=20; %minimum pressure

%% And do the plotting

%Loop through all files and profiles; plot each with an offset
cc=0; %profile counter
figure(1); clf; hold on; grid on; axis ij; box on; zoom on; ax=gca;
for ii=1:length(fnames);
    %Load file
    fprintf('Loading %s...\n',fnames{ii});
    load([pth fnames{ii}]);
    disp('Done loading. Plotting...');
    %Loop through each profile in file and plot
    for jj=1:length(Microrider)
        cc=cc+1; %Keep track of profiles
        %First the glider data
        iplot = Glider(jj).P>minP; %only where P>minP
        plot(Glider(jj).CT(iplot)+cc*delT, Glider(jj).P(iplot), 'k:');
        %Then the microrider data
        micP = Microrider(jj).P(1:scDn:end); %scaled down resolution
        micT1 = Microrider(jj).micT1(1:scDn:end); %scaled down resolution
        micT2 = Microrider(jj).micT2(1:scDn:end); %scaled down resolution
        iplot = micP>=minP; %only where P>minP
        plot(micT1(iplot)+cc*delT, micP(iplot), 'linewidth', 1.5)
        plot(micT2(iplot)+cc*delT, micP(iplot), 'linewidth', 1.5)
        ax.ColorOrderIndex=1; %reset colour order for plotting
        [~, itxt] = max(micP);
        %add absolute cast number
        text(micT1(itxt)+cc*delT, micP(itxt)+5, num2str(cc), 'fontsize', 12)
        %add cast number of this file
        text(micT1(itxt)+cc*delT, micP(itxt)+15, sprintf('%d-%d',ii,jj),...
            'fontsize', 12)
    end
    drawnow; disp('Done plotting');
end
clear micP micT1 micT2 iplot ii jj ax itxt





































