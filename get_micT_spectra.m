%parent: t02_make_noise_model.m

function output = get_micT_spectra(data,params,kmax)
%
%Calculate and integrate the temperature gradient spectra so we can develop
%a noise model for the MicT probes. This is just a very paired down version
%of micT_eps_chi.m
%
%Dependencies:
%   :closeto()
%
%B.Scheifele 2017-08

%Constants
q = 3.4; %Batchelor spectrum parameter
Dt = 1.44e-7; %molecular temperature diffusivity [m^2/s]
maxFreq = 98; %Use anti-aliasing frequency for highest frequency. RSI TN-10 

%Parse spectral input parameters
overlap = params.overlap; %overlap between consecutive eps/chi estimates
fftLen = params.fftLen; %number of samples per FFT
dissLen = params.dissLen; %number of consecutive samples per eps/chi estimate
freqInstr = params.freq; %instrument frequency, exact, for FFT

%Calculate the number of dissipation estimates that can be made
nEstimates = 1 + floor((length(data.mtime)-dissLen) / (dissLen-overlap));
%Spacing and size of the output frequency vector
dFreq = freqInstr/fftLen; %spacing of frequency vector
nFr = floor(fftLen/2); %size of full frequency vector
nFrSt = floor(maxFreq/dFreq); %size of shortened frequency vector

%Pre-allocate a few matrices. Two micT channels
output.mtime = nan(nEstimates,1); %matlab serial time
output.P = nan(nEstimates,1); 
output.U = nan(nEstimates,1);
output.f = nan(nEstimates,nFrSt); %frequency [Hz = 1/s]
output.spec_f = nan(nEstimates,2,nFrSt); %spectrum in f space [K^2 s]
output.k = nan(nEstimates,nFrSt); %wavenumber [cpm = 1/m]
output.spec_k = nan(nEstimates,2,nFrSt); %gradT spectrum, k space [K^2/cpm]
output.variance = nan(nEstimates,2);

%Loop through Chi estimates
for ii=1:nEstimates
    
    %indices for this estimate
    inds = [1:dissLen] + (dissLen-overlap)*(ii-1);
    
    %Data for this estimate
    T = nan(length(inds),2);
    T(:,1) = data.micT1(inds);
    T(:,2) = data.micT2(inds);
    P = data.P(inds);
    U = data.U(inds);
    meanTime = mean(data.mtime(inds),'omitnan');
    meanU = mean(U,'omitnan');
    meanP = mean(P,'omitnan');
       
    %Skip if there is no valid U measurement
    if ~isfinite(meanU)
        fprintf(2,'skipping datum %d, no U available\n',ii);
        continue  %if no U, return to top of loop and iterate
    end
    
    %Some output variables
    output.mtime(ii) = meanTime;
    output.P(ii) = meanP;
    output.U(ii) = meanU;

    %calculate spectra using ODAS. 50% overlap. Cos Window. Linear detrend.
    %<T^2> = int{S(f)}df, so units of S(f) are [K^2 s]
    S_fr = nan(2,nFr+1); fr=nan(1,nFr+1);
    for jj=1:2
        [S_fr(jj,:),fr(:)] = csd_matrix_odas(T(:,jj), [], fftLen,...
            freqInstr, [], fftLen/2, 'linear');
    end
    %First datapoint is always (0,0). Don't want this...
    S_fr = S_fr(:,2:end); fr = fr(2:end);
        
    %Remove highest frequencies, where anti-aliasing filter kicks in
    ikeep = fr<maxFreq; fr = fr(ikeep); S_fr = S_fr(:,ikeep);
    
    %Don't correct for slow sensor response here - note, this is
    %different than for the actual calculations in the next step

    %Keep these
    output.f(ii,:) = fr;
    output.spec_f(ii,:,:) = S_fr;
    
    %Calculate wavenumber vector. Assume Taylor's frozen turbulence hypoth
    k=fr/meanU; %[cpm] = [1/m]
    
    %From S(f) to S(k). Multiply by U to preserve the total variance.
    %You can see this from the change of variables: if f=kU then df=Udk and
    %and int{S(f)}df = int{S(k)U}dk. Assume frozen turbulence.
    S_k = S_fr * meanU;
    
    %Calculate the gradient power spectrum. If <T^2> = int{S(k)}dk then
    %<(dT/dz)^2> = int{k^2 S(k)}dk. Units are [K^2 cpm] or [K^2/m]. Note
    %that our k is in cpm so to get the derivative right, we need to
    %include the factor of 2pi: If d/dx(sin(kx)) = k cos(kx) then
    %that k must be in rad/m
    S_k = (2*pi*[k;k]).^2 .* S_k;
    
    %Keep these
    output.k(ii,:) = k;
    output.spec_k(ii,:,:) = S_fr;
    
    %Loop through two probes
    for jj=1:2
        
        %Get Spectra for this probe
        Sobs = S_k(jj,:); %observed
        
        %Back to top of loop if there was problem calculating spectrum
        if any(~isfinite(Sobs))
            fprintf(2,'skipping datum %d, spectrum not strictly finite\n',ii);
            output.spec_f(ii,jj,:)=nan;
            output.spec_k(ii,jj,:)=nan;
            continue
        end
        
        %Indices over which to integrate
        kInds = 1:closeto(k,kmax);
        
        %Integrate
        variance = trapz(k(kInds),Sobs(kInds));
        
        %remove if integral is nonfinite
        if ~isfinite(variance) 
            fprintf(2,'skipping datum %d, integral not strictly finite\n',ii);
            variance = nan;
            output.spec_f(ii,jj,:)=nan;
            output.spec_k(ii,jj,:)=nan;
            continue 
        end
        
        output.variance(ii,jj)=variance;

    end %end loop over probes
    
end %end loop over estimates

end %end function






