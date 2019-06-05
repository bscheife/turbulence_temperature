
%Here we calculate temperature gradient spectra from all available data,
%and then average the lowest 1% of spectra to create our estimate for the
%noise spectrum. We need the noise spectrum for the MLE approach (in the
%next processing step) to get the Batchelor wavenumber. See Ruddick 2000.

%Runtime is ~20 mins for Geotraces2015 data

%The spectral processing parameters must match those you're going to use
%for the actual calculations in the next step

clear
tic

%% %%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%

pth = '../02_collect/'; %path to input file(s)
fnames = {'GlidMic01.mat'; 'GlidMic02.mat'; 'GlidMic03.mat'}; %input files
% fnames = {'GlidMicTest.mat'};
out_fname = 'noise_spectra_geotr2015.mat'; %output filename for noise spect
save2mat = false;

%Spectral processing parameters
fftLen = 4; %Seconds, length of one fourier transform. Integer [4]
nFFT = 10; %Number of FFT lengths to use for each dissipation estimate. Integer
freq = 512; %Frequency of shear record (Hz, integer). 512 for Microrider
%Fer et al '14 and Peterson & Fer 2014 use nSecs=4, nFFT=3;

%maximum k upto which to integrate. We can just use a constant here...
kmax = 30; %[cpm]

%% %%%%%%%%%%%%% COUNT & INITIALIZE %%%%%%%%%%%%%%%%%%

%Spectral processing parameters
params.fftLen = fftLen*freq; %FFT length, number of samples
params.dissLen = nFFT*params.fftLen; %number of samples in epsilon estimate
params.overlap = params.dissLen/2; %overlap between consec chi calcs: 50%

%count total number of casts - yes, it's faster than not pre-allocating!
disp('Counting all casts...')
nCasts=0;
for ii=1:length(fnames)
    load([pth fnames{ii}],'Glider');
    nCasts=nCasts+length(Glider);
end
disp(['Done counting: ' num2str(nCasts)]);

%initialize output variable
spectra(nCasts).k=[];

%clean up
clear ii fftLen nFFT freq

%% %%%%%%%%% BEGIN NESTED LOOP W/ CALCULATION %%%%%%%%%%%%%

cc=0; %cast counter
%loop over all files. They're big, so only load one at a time

for ii=1:length(fnames)    
    
    %Load data file
    disp(['Loading file: ' fnames{ii}]);
    load([pth fnames{ii}]);
    disp('Done loading. Processing...');
    
    %loop over casts in this file
    for jj=1:length(Microrider)        
    
        %Count
        cc=cc+1;
        disp(['Cast Number ' num2str(cc) '...']);
        
        %Measurements needed for spectral calculations
        params.freq = Microrider(jj).frequ; %exact measurement frequency
        data.micT1 = Microrider(jj).micT1; %micT record 1
        data.micT2 = Microrider(jj).micT2; %micT record 2
        data.P = Microrider(jj).P; %pressure record
        data.U = Microrider(jj).U; %microrider speed through water
        data.mtime = Microrider(jj).mtime; %matlab time   
        
        %Skip calculations if no MicT record, return to top of loop
        if isempty(data.micT1) & isempty(data.micT2)
            warning('No micT measurements in this cast. Skipping');
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function to Calculate Spectra and Chi
        output = get_micT_spectra(data,params,kmax);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Collect all values into our structure array
        fldnames = fieldnames(output);
        for kk=1:length(fldnames)
            spectra(cc).(fldnames{kk}) = output.(fldnames{kk});
        end
                
    end %end loop over casts.        
    clear casts jj
    
end %end loop over files
clear ii cc
toc

%%

%all Variance in a 2d array
allVar = cell2mat({spectra.variance}');

%sum two Variance estimates & find indeces of lowest 1% of chi estimates.
%Summing ensures simultaneous estimates are both low and finite
N = floor(0.01*length(allVar));
disp(['N is ' num2str(N)]);
[SortSumVar idSort]=sort(sum(allVar,2),'missingplacement','last');
id = idSort(1:N); %indeces of the lowest 1 percent
%id = idSort(50*N:51*N); %indeces of the central 1 percent

clear SortSumChi idSort

%put all spectra in a single 3D array
allSpec_f = cell2mat({spectra.spec_f}');
%frequency vector is the same for all of them
fr = unique(cell2mat({spectra.f}'),'rows');
if ~isrow(fr)
    error('looks like the errors arent unique')
end
%now isolate the lowest N spectra
spec_f = allSpec_f(id,:,:);

%Plot the power spectra in fr space
figure(1); clf; cm = get(groot,'DefaultAxesColorOrder');
hold on; box on; grid on; zoom on;
set(gca,'xscale','log','yscale','log');
for ii=1:2
    for jj=1:N
        ph=plot(fr,squeeze(spec_f(jj,ii,:))',...
            'linewidth',2,'color',cm(ii,:));
        ph.Color(4)=0.075;
    end
end
xlabel('fr (1/s)'); ylabel('S(fr)')

%This looks promising. The two probes differ substantially, indicating we
%really are in the noise since they appear to have different noise floors

%The frequency vectors are all on the same grid, so we can just average!
%But we should do it in log space
S_noise = squeeze(10.^mean(log10(spec_f)));

%And plot
plot(fr,S_noise,'color','k','linewidth',2)

%% And save

fr_noise = fr;
if save2mat
    save(out_fname,'fr_noise','S_noise');
end


































