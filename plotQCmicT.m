%parent: chooseQCmicT.m and t05_micT_QC.m

function plotQCmicT(MicTSpectra, fldname, datalims);
%make figure for chooseQCmicT.m. third argument (datalims) is optional.

if nargin<2, error('too few inputs'); end

%x and y axes
n = 1:length(MicTSpectra);
pmax = max(cell2mat({MicTSpectra.P}'));
p = [0:0.1:pmax]; 

%Make matrices
[dataMat1 dataMat2] = Fields2Matrix(MicTSpectra, fldname, 'P', p);
dataCell = {dataMat1 dataMat2};

%if no input for data limits
if nargin<3
    mindat = min(log10([dataMat1(:); dataMat2(:)]));
    maxdat = max(log10([dataMat1(:); dataMat2(:)]));
    datalims = [mindat maxdat];
end

for jj=1:2
    dataMat = log10(dataCell{jj});
    %Section, top left
    ax1=subtightplot(8,4,(jj-1)*16+[1:3 5:7 9:11],[],[0.08 0.05],0.05); 
    ax1.XTickLabel=[];
    box on; grid on; axis ij; hold on; zoom on;
    pcolor(n,p,dataMat); shading flat; colorbar;
    ylim([0 pmax+5]); xlim([min(n) max(n)]);
    set(gca,'color',0*[1 1 1]); colorbar; caxis(datalims);
    if jj==1, title(['log10 ' fldname]), end; 
    ylabel('p dbar')
    %Mean profile, top right
    ax2=subtightplot(8,4,(jj-1)*16+[4:4:12],[],[0.08 0.05],0.05);
    set(ax2,'xticklabel',[],'yticklabel',[]);
    box on; grid on; hold on; zoom on; axis ij;
    plot(mean(dataMat,2,'omitnan'),p,'.','linewidth',2);
    ylim([0 pmax+5]); xlim(datalims)
    if jj==1, title(['log10 ' fldname]), end;
    %Mean time series
    ax3=subtightplot(8,4,(jj-1)*16+[13:15],[],[0.08 0.05],0.05);
    if jj==1, set(ax3,'xticklabel',[]), end;
    box on; grid on; hold on; zoom on;
    plot(n, mean(dataMat,'omitnan'),'.', 'linewidth', 2)
    xlim([min(n) max(n)]); ylim(datalims)
    ax3.Position(3)=ax1.Position(3);
    ylabel(fldname);
    if jj==2, xlabel('cast number'); end
    %Histogram
    ax4=subtightplot(8,4,(jj-1)*16+16,[],[0.08 0.05],0.05); 
    box on; grid on; hold on; zoom on;
    ax4.YTickLabel=[];
    if jj==1, ax4.XTickLabel=[];, end
    histogram(dataMat,25);
    xlim(datalims)
    if jj==2, xlabel(fldname), end;
end




end