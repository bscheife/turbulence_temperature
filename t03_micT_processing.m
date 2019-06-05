

clear
% clearvars -except Glider Microrider
close all
tic

%% %%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%

%Input data file(s). Specify filenames as: {'fname1.mat'; 'fname2.mat'}
in.fnames = {'GlidMic01.mat'; 'GlidMic02.mat'; 'GlidMic03.mat'}; %file name(s)
% in.fnames = {'GlidMicTest.mat'};
in.pth = '../02_collect/'; %path to prepend to input file(s)

%Debug options. Useful if looking at specific casts and loading no more
%than one file at a time (i.e. if length(in.fnames)==1). If batch processing,
%use [recommended] values
in.preLoaded = false; %true/[false]. If true, assumes file already loaded. Quicker.
in.whichCasts = 'all'; %string: individual ('5'), range ('1:10'), or [all] ('all')

%Output options
out.saveFile = true; %[true]/false, save an output MAT file.
% out.fname = 'MicT_preQC.mat'; %output filename for combined micT epsilon
% out.fname = 'MicT_preQC_qB2.mat'; %to test sensitivity to qB
out.fname = 'MicT_preQC02.mat'; %with new Chi definition 2018-07-26

%Draw diagnostic figures?
out.drawFigs = false; %true/[false]. Use false when batch processing
out.nPause = -10; %stop every Nth chi estimate to draw spectra. negative for none

%Spectral processing parameters
fftLen = 4; %Seconds, length of one fourier transform. Integer [4]
nFFT = 10; %Number of FFT lengths to use for each dissipation estimate. Integer
freq = 512; %Frequency of shear record (Hz, integer). 512 for Microrider
%Fer et al '14 and Peterson & Fer 2014 use nSecs=4, nFFT=3;

%% %%%%%%%%%%%%% COUNT & INITIALIZE %%%%%%%%%%%%%%%%%%

%Some universal spectral processing parameters. 
spectr.fftLen = fftLen*freq; %FFT length, number of samples
spectr.dissLen = nFFT*spectr.fftLen; %number of samples in epsilon estimate
spectr.overlap = spectr.dissLen/2; %overlap between consec chi calcs: 50%

%count total number of casts - yes, it's faster than not pre-allocating!
disp('Counting all casts...')
nCasts=0;
for ii=1:length(in.fnames)
    load([in.pth in.fnames{ii}],'Glider');
    nCasts=nCasts+length(Glider);
end
disp(['Done counting: ' num2str(nCasts)]);

%initialize output variable
MicTSpec(nCasts).mtime=[];

clear ii fftLen nFFT freq

%% %%%%%%%%%%%%%% BEGIN NESTED LOOP %%%%%%%%%%%%%%%%%%

cc=0; %cast counter
%loop over all files. They're big, so only load one at a time
for ii=1:length(in.fnames)    

    %Load data file
    if ~in.preLoaded
        clear Glider Microrider
        disp(['Loading file: ' in.fnames{ii} '. Might take a few minutes...']);
        load([in.pth in.fnames{ii}]);
        disp('Done loading. Processing...');
    end
    
    %Parse which casts are requested
    if strcmpi(in.whichCasts,'all')
        casts=1:length(Microrider);
    else
        casts=eval(in.whichCasts);
    end
    
    %loop over casts in this file
    for jj=casts        
    
        %Count
        cc=cc+1;
        disp(['Cast Number ' num2str(cc) '...']);
        
        %Measurements needed for spectral calculations
        spectr.freq = Microrider(jj).frequ; %exact measurement frequency
        data.micT1 = Microrider(jj).micT1; %micT record 1
        data.micT2 = Microrider(jj).micT2; %micT record 2
        data.P = Microrider(jj).P; %pressure record
        data.U = Microrider(jj).U; %microrider speed through water
        data.mtime = Microrider(jj).mtime; %matlab time   
        data.pitch = Microrider(jj).pitch; %pitch
        data.aoa = Microrider(jj).aoa; %angle of attack
        
        %Skip calculations if no MicT record, return to top of loop
        if isempty(data.micT1) & isempty(data.micT2)
            warning('No micT measurements in this cast. Skipping');
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function to Calculate Spectra, Chi, and Epsilon
        EpsChi = micT_eps_chi(data,spectr,out);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Collect all values into our structure array
        fldnames = fieldnames(EpsChi);
        for kk=1:length(fldnames)
            MicTSpectra(cc).(fldnames{kk}) = EpsChi.(fldnames{kk});
        end
                
    end %end loop over casts.        
    clear casts jj
    
end %end loop over files
clear ii cc
toc

%Save if requested
if out.saveFile
    tic
    save(out.fname,'MicTSpectra');
    toc
end






































