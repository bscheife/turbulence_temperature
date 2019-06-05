%parent: t05_micT_QC.m

function imageQCmicT(MicTSpectra, fldname, orig_fldname, datalims);
%make figure to display what this QC condition does to the data. 
%fourth argument (datalims) is optional.

if nargin<3, error('too few inputs'); end

%x and y axes
n = 1:length(MicTSpectra);
pmax = max(cell2mat({MicTSpectra.P}'));
p = [0:0.1:pmax]; 

%Make matrices
[dataOrig1 dataOrig2] = Fields2Matrix(MicTSpectra, orig_fldname, 'P', p);
[dataQC1 dataQC2] = Fields2Matrix(MicTSpectra, fldname, 'P', p);
dataOrig_cell = {dataOrig1 dataOrig2};
dataQC_Cell = {dataQC1 dataQC2};

%if no input for data limits
if nargin<4
    mindat = min(log10([dataQC1(:); dataQC2(:)]));
    maxdat = max(log10([dataQC1(:); dataQC2(:)]));
    datalims = [mindat maxdat];
end

axcount = 0;
for jj=1:2
    %get data
    dataQC = log10(dataQC_Cell{jj});
    dataOrig = log10(dataOrig_cell{jj});
    %original left
    axcount = axcount+1;
    ax(axcount)=subtightplot(2,2,(jj-1)*2+[1],[],[0.08 0.05],0.05);
    box on; grid on; axis ij; hold on; zoom on;
    pcolor(n,p,dataOrig); shading flat; colorbar; 
    set(gca,'color','k'); caxis(datalims)
    if jj==1, 
        title('Original');
        set(gca,'xticklabel', []);
    else
        xlabel('cast number');
    end;
    ylim([0 pmax+5]); xlim([min(n) max(n)]);

    %QC'd right
    axcount = axcount+1;
    ax(axcount)=subtightplot(2,2,(jj-1)*2+[2],[],[0.08 0.05],0.05);
    box on; grid on; axis ij; hold on; zoom on;
    pcolor(n,p,dataQC); shading flat; colorbar; 
    set(gca,'color','k'); caxis(datalims)
    set(gca,'yticklabel',[]);
    if jj==1 
        title(fldname);
        set(gca,'xticklabel',[]);
    else
        xlabel('cast number');
    end
    ylim([0 pmax+5]); xlim([min(n) max(n)]);

end

% ax(3).Position(3)=ax(1).Position(3);
% ax(4).Position(3)=ax(2).Position(3);
















end