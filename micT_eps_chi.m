
function epschi = micT_eps_chi(data,spectr,out)
%Use Ruddick's method for calculating chi and epsilon from micT. The
%overall algorithm is described well by Steinbuck et al 2009. The MLE for
%estimating kB (and therefore epsilon) is described by Ruddick et al 2000. 
%
%with codes provided by J.Carpenter and B.Ruddick
%
%Dependencies:
%   :closeto()
%   :micT_noise()
%   :fit_kB()
%       :Batch_spec()
%       :C11()
%       :integrate_Batchelor()
%
%B.Scheifele 2017-05

%2018-07-26 B.Sch. Make the calculation of chi less restrictive.
%Previously, we used X_T as our best estimate of chi. Now we're using X_T
%only for the epsilon calculation, and using the new variable X_int for the
%best estimate of chi. The advantage is that we don't need to make
%assumptions about Batchelor-like shapes to get an estimate of chi - we
%just integrate the observed spectrum. Note, you can still get X_T in the
%output if you want it: just add X_lw + X_obs + X_hw.


%Constants
q = 3.4; %Batchelor spectrum parameter % q = 3.7; %alternate qB value
Dt = 1.44e-7; %molecular temperature diffusivity [m^2/s]
maxFreq = 98; %Use anti-aliasing frequency for highest frequency. RSI TN-10 

%Parse spectral input parameters
overlap = spectr.overlap; %overlap between consecutive eps/chi estimates
fftLen = spectr.fftLen; %number of samples per FFT
dissLen = spectr.dissLen; %number consecutive samples per eps/chi estimate
freqInstr = spectr.freq; %instrument frequency, exact, for FFT

%degrees of freedom is 2xnSpectra, where nSpectra is the number of
%spectra averaged to get one "useable" spectrum. This may be too naive
%but I don't know a better way to get the DOF.
dof = 2 * (2*dissLen/fftLen-1);

%Parse figure options
drawFigs = out.drawFigs; %draw figures?
nPause = out.nPause; %pause how often?
cm=get(groot,'DefaultAxesColorOrder'); %default color map
global stop2draw; %global variable for choosing whether to draw figs

%Calculate the number of dissipation estimates that can be made
nEstimates = 1 + floor((length(data.mtime)-dissLen) / (dissLen-overlap));
%Spacing and size of the output frequency vector
dFreq = freqInstr/fftLen; %spacing of frequency vector
nFr = floor(fftLen/2); %size of full frequency vector
nFrSt = floor(maxFreq/dFreq); %size of shortened frequency vector

%Pre-allocate a few matrices. For now I'll just assume two micT channels
epschi.mtime = nan(nEstimates,1); %matlab serial time
epschi.P = nan(nEstimates,1); 
epschi.U = nan(nEstimates,1);
epschi.nu = nan(nEstimates,1);
epschi.dTdz = nan(nEstimates,1); %mean vertical gradient over chi estimate
epschi.dTdx = nan(nEstimates,1); %mean along path gradient ovr chi estimate
epschi.eps = nan(nEstimates,2);
epschi.chi = nan(nEstimates,2);
epschi.cox = nan(nEstimates,2); %cox number
epschi.k = nan(nEstimates,nFrSt); %wavenumber [cpm = 1/m]
epschi.spec_measu = nan(nEstimates,2,nFrSt); %measured spectrum [K^2/cpm]
epschi.spec_noise = nan(nEstimates,2,nFrSt); %noise model spectrum [K^2/cpm]
epschi.spec_batch = nan(nEstimates,2,nFrSt); %batchelor spectrum [K^2/cpm]
epschi.spec_Y_ratio = nan(nEstimates,2,nFrSt); %ratio S_obs/(S_Batch+noise)
epschi.kmin = nan(nEstimates,2); %lowest k used for integral and Batch fit
epschi.kmax = nan(nEstimates,2); %largest k used for integral and Batch fit
epschi.kB = nan(nEstimates,2); %Batchelor wavenumber
epschi.kstar = nan(nEstimates,2); %used for choosing kmin
epschi.num_k = nan(nEstimates,2); % #-of-k's used for integrl and Batch fit
epschi.X_obs = nan(nEstimates,2); %directly observed portion of X
epschi.X_lw = nan(nEstimates,2); %lower wavenumber correction to X
epschi.X_hw = nan(nEstimates,2); %upper wavenumber correction to X
epschi.MAD = nan(nEstimates,2); %mean absolute deviation observed to theor
epschi.dof = dof; %degrees of freedom for spectrum = 2x # of periodograms

%The Cox number (epschi.cox) is a physical quantity, the ratio of the
%variance in the turbulent field to that in the mean field:
% 3<(dT'/dx)^2> / (dT/dx)^2
%Since we detrend before calculating (dT'/dx)^2, Cox approaches zero for
%laminar flow. If we did not detrend, it should approach one for laminar
%flow. Unlike for the \chi calculation, we want to use all the available
%variance and integrate from the lowest available wavenumber to the noise

%Loop through Chi and Epsilon estimates
for ii=1:nEstimates
    
    %disp(ii);
    
    %Parse flag to draw figures
    if drawFigs && mod(ii,nPause)==0 && nPause>0
        stop2draw=true;
    else
        stop2draw=false;
    end
    
    %indices for this estimate
    inds = [1:dissLen] + (dissLen-overlap)*(ii-1);
    
    %Data for this estimate
    T = nan(length(inds),2);
    T(:,1) = data.micT1(inds);
    T(:,2) = data.micT2(inds);
    P = data.P(inds);
    U = data.U(inds);
    time = datetime(datevec(data.mtime(inds))); %for plotting only
    meanTime = mean(data.mtime(inds),'omitnan'); %for output
    meanU = mean(U,'omitnan');
    meanP = mean(P,'omitnan');
    pitch = mean(data.pitch(inds),'omitnan'); %glider pitch, ccw positive
    aoa = mean(data.aoa(inds),'omitnan'); %glider AOA, ccw positive
    gamma = pitch - aoa; %glide angle, |gamma|>|pitch| always
    mean_dTdz = mean((T(end,:)-T(1,:)) ./ (P(end)-P(1))); %vertical grad
    mean_dTdx = mean_dTdz * sind(-gamma); %along-path temperature gradient
    nu = visc35(mean(T(:),'omitnan')); %kinematic viscosity m2/s
       
    %Skip if there is no valid U measurement
    if ~isfinite(meanU)
        fprintf(2,'skipping datum %d, no U available\n',ii);
        continue  %if no U, return to top of loop and iterate
    end
    
    %Some output variables
    epschi.mtime(ii) = meanTime;
    epschi.P(ii) = meanP;
    epschi.U(ii) = meanU;
    epschi.nu(ii) = nu;
    epschi.dTdz(ii) = mean_dTdz;
    epschi.dTdx(ii) = mean_dTdx;
        
    %calculate spectra using ODAS. 50% overlap. Cos Window. Linear detrend.
    %<T^2> = int{S(f)}df, so units of S(f) are [K^2 s]
    S_obs = nan(2,nFr+1); fr=nan(1,nFr+1);
    for jj=1:2
        [S_obs(jj,:),fr(:)] = csd_matrix_odas(T(:,jj), [], fftLen,...
            freqInstr, [], fftLen/2, 'linear');
    end
    %First datapoint is always (0,0). Don't want this...
    S_obs = S_obs(:,2:end); fr = fr(2:end);
    
    %Remove highest frequencies, where anti-aliasing filter kicks in
    ikeep = fr<maxFreq; fr = fr(ikeep); S_obs = S_obs(:,ikeep);
    
    %Plotting
    if stop2draw
        figure(1); clf;
        %Plot temperature as function of time
        ax(1)=subplot(241); grid on; box on; hold on; axis ij; zoom on;
        plot(T,time); xlabel('T'); ylabel('time');
        %Temperature variance power spectrum, frequency space
        ax(2)=subplot(2,4,[2:4]); grid on; box on; hold on; zoom on;
        set(gca,'yscale','log','xscale','log','yaxislocation','right');
        xlabel('fr (Hz)'); ylabel('S_T  (K^2/Hz)'); title('Power Spectrum')
        plot(fr, S_obs);
        %Temperature as function of pressure
        ax(3)=subplot(245); grid on; box on; hold on; axis ij; zoom on;
        plot(T,P); xlabel('T'); ylabel('P (dbar)')
    end
    
    %Correct for slow sensor response at high frequencies with transfer
    %function H^2. We use Sommer et al 2013 which is independent of U, but
    %could also use Gregg & Meager where tau=0.007*U^(-0.32)
    tau = 0.010; %[s]
    Hsq=((1+(2*pi*[fr;fr]*tau).^2).^2).^(-1); %transfer fn
    S_obs=S_obs./Hsq; %corrected, [K^2 s]

    %Noise spectrum in frequency space, corrected for slow response
    S_ns=micT_noise(fr)./Hsq; %units [K^2 s]
    
    if stop2draw
        figure(1); axes(ax(2)); ax(2).ColorOrderIndex=1;
        %Plot corrected noise spectrum in frequency space
        plot(fr, S_ns,'--'); ax(2).ColorOrderIndex=1;
        %Plot corrected temperature variance power spectrum in frq space
        plot(fr, S_obs, 'linewidth',3);
    end
    
    %Calculate wavenumber vector. Assume Taylor's frozen turbulence hypoth
    k=fr/meanU; %[cpm]  (has dimension 1/L)
   
    %From S(f) to S(k). Multiply by U to preserve the total variance.
    %You can see this from the change of variables: if f=kU then df=Udk and
    %and int{S(f)}df = int{S(k)U}dk. Assume frozen turbulence. Units are
    %[K^2/cpm]
    S_obs = S_obs * meanU;
    S_ns = S_ns * meanU;
    
    %Calculate the gradient power spectrum. If <T^2> = int{S(k)}dk then
    %<(dT/dz)^2> = int{k^2 S(k)}dk (if k is in units of rad/m). Units are
    %[K^2 cpm]. Note that our k is in cpm so to get the derivative right,
    %we need to include the factor of 2pi: If d/dx(sin(kx)) = k cos(kx)
    %then that k needs to have units rad/m for this to work, but our k has
    %units cpm, so we need to multiply by 2pi rad/cyc. Stupid f*** 2pi
    %factor again
    %Update - checked w/ Rolf, this is right (04/2018)
    S_obs = (2*pi*[k;k]).^2 .* S_obs;
    S_ns = (2*pi*[k;k]).^2 .* S_ns;

    %Some output variables
    epschi.k(ii,:)=k;
    epschi.spec_measu(ii,:,:) = S_obs; %measured spectrum [K^2/cpm]
    epschi.spec_noise(ii,:,:) = S_ns; %noise model spectrum [K^2/cpm]
        
    %Plotting gradient power spectrum in k space
    if stop2draw
        figure(1);
        %Temperature gradient power spectrum, as function of k
        ax(4)=subplot(2,4,6:8); hold on; grid on; box on; zoom on;
        set(gca,'yscale','log','xscale','log','yaxislocation','right');
        ylabel('S  (K^2/cpm)'); xlabel('k (cpm)'); title('Gradient Power Spectrum')
        lh=plot(k,S_obs,'linewidth',3); 
        ax(4).ColorOrderIndex=1;
        %Corresponding noise spectrum
        plot(k, S_ns, '--');
    end
    
    %Loop through two probes
    for jj=1:2
        
        %Get Spectra for this probe
        Sobs = S_obs(jj,:); %observed
        Snos = S_ns(jj,:); %noise     
        %Back to top of loop if there was problem calculating spectrum
        if any(~isfinite(Sobs))
            fprintf(2,'skipping datum %d, spectrum not strictly finite\n',ii);
            continue
        end
        
        %Choose lower and upper integration limits to get Chi
        %Lower limit is just an inital guess. This gets revisited later
        kL_ind=1;
        %Upper limit is intersection of measured and 2x noise spectrum
        ind_U=find(Sobs<2*Snos,1);
        if isempty(ind_U)
            kU=k(end);
            kU_ind=length(k);
        else
            kU_ind=ind_U;
            kU=k(kU_ind);
        end
        clear indU
        %Indices over which to integrate
        kInds = kL_ind:kU_ind;
        
        %%%%%%%% Interlude - get directly-integrated Chi %%%%%%%%%
        %Calculate X_int. First find upper integration limit - constrain to
        %be at least 7 cpm. Note, we have to include the zero point here
        %when integrating because, unlike for the later chi calculation,
        %there's no correction term for low wavenumbers (i.e. no X_lw).
        if kU<7
            iupper = find(k>7,1);
            fprintf(2,'using kmax = 7. Datum %d\n', ii);
        else
            iupper = kU_ind;
        end
        ks = [0 k(1:iupper)]; %include zero here
        Sp = Sobs(1:iupper)-Snos(1:iupper); 
        Sp = [0 Sp]; %include zero here
        X_int = 6*Dt*trapz(ks,Sp);
        %if this is less than zero, just integrate the noise spectrum
        if X_int<0
            X_int = 6*Dt*trapz(ks,[0 Snos(1:iupper)]);
        end
        %clean up
        clear ks Sp
        %Get the cox number
        Cox = (X_int/2/Dt) / (mean_dTdx).^2;
        if isfinite(X_int)
            epschi.chi(ii,jj) = X_int; %keep this value
            epschi.cox(ii,jj) = Cox; %keep this value
        else
            disp('Couldnt get chi by direct integration. Skipping')
            continue %back to top of loop, skip the rest
        end
        %%%%%%%%%%%%%%%%%%% End Interlude %%%%%%%%%%%%%%%%%%
        
        if kU_ind==1
            %If the upper index is the first k, we're too close to the 
            %noise to continue. Let's turn around
            fprintf(2,'skipping datum %d, unable to get upper integration limit\n', ii);
            continue
        end
        
        %Integrate (measured)-(noise) for initial estimate of X_T.
        %Factor of 6 from assuming isotropy
        X_T = 6*Dt*trapz(k(kInds),Sobs(kInds)-Snos(kInds));
        
        %top of loop if X_T (i.e. the integral) is non-finite
        if ~isfinite(X_T) 
            continue %back to top of loop
        end
        
        %Try to calculate best fit Batchelor Spectrum using MLE, following
        %Ruddick et al 2000. Exit if this fails
        try
            [~,kB,f11]=fit_kB(k(kInds), Sobs(kInds), Snos(kInds),...
                                       X_T, nu, Dt, dof, q);
        catch
            fprintf(2,'skipping datum %d, unable to get batchelor fit\n', ii)
            continue %back to top of loop
        end
        
        %calculate epsilon from batchelor spectrum. Theoretical expression
        %requires kB to have units of rad/m. Eq.3 in Ruddick 2000
        epsilon=(kB*2*pi)^4*nu*Dt^2; 
        
        %Iterate twice, modifying kL to avoid contamination from fine
        %structure etc. Use kB to calculate k_star (upper limit of inertial
        %convective subrange) and then choose the larger of k_star and k(1)
        itrErrorFlag=false;
        for iteration=1:2
            %Define k_star as in Steinbuck et al 2009
            k_star=0.04*(nu/Dt)^(-1/2)*kB;
            %For kL, choose the larger of k(1) and 3*k_star. I found this
            %gives better Batchelor fits than comparing k(1) and 1*k_star
            if 3*k_star>=k(1)
                %be conservative and round up
                kL_ind=find(k>3*k_star,1);
                kL=k(kL_ind);
            else
                kL_ind=1;
                kL=k(kL_ind);
            end
            %Indices over which to integrate
            kInds = kL_ind:kU_ind;
            %Exit if there are fewer than 3 useable points
            if length(kInds)<3
                fprintf(2,'skipping datum %d: fewer than 3 useable points\n',ii);
                itrErrorFlag=true;
                continue
            end
            %Integrate using new limits
            X_obs = 6*Dt*trapz(k(kInds),Sobs(kInds)-Snos(kInds));
            %Find correction terms for unresolved portions of spectrum
            X_lw = 6*Dt*integrate_Batchelor(0,kL,kB,X_T,Dt,q);
            X_hw = 6*Dt*integrate_Batchelor(kU,100*kU,kB,X_T,Dt,q);
            X_T = X_lw+X_obs+X_hw; %total new X
            %Calculate best-fit Batchelor Spectrum using MLE
            try
                [~,kB,f11]=fit_kB(k(kInds), Sobs(kInds), Snos(kInds),...
                    X_T, nu, Dt, dof, q);
            catch
                fprintf(2,'skipping datum %d: unable to get batchelor fit\n',ii);
                itrErrorFlag=true;
                continue
            end
            %Calculate the full Batchelor spectrum, not just fitted part
            f11 = Batch_spec(kB,k,X_T,Dt,q);
            %Calculate refined epsilon
            epsilon=(kB*2*pi)^4*nu*Dt^2;
        end
        
        %If there was a problem in the previous loop, exit
        if itrErrorFlag
            continue
        end
        
        if stop2draw
            %output kstar to the screen
            fprintf('k_star is %.2f cpm\n',k_star);
            %plot batchelor spectrum and lower and upper limits used for
            %integration and batchelor fitting
            miny=min(S_ns(:)); maxy=max(S_ns(:));
            minx=0.9*min(k); maxx=max(k);
            linestyle={':', '--'};
            figure(1); axes(ax(4)); ylim([miny maxy]); xlim([minx maxx]);
            plot(k, f11, 'color', cm(jj,:))
            plot(kL*[1 1], [miny maxy], linestyle{jj}, 'color', cm(jj,:),'linewi',2)
            plot(kU*[1 1], [miny maxy], linestyle{jj}, 'color', cm(jj,:),'linewi',2)
            uistack(flipud(lh),'top'); clear miny maxy minx maxx
        end
        
        %ratio observed to theoretical
        Y = Sobs(kInds)./(f11(kInds)+Snos(kInds)); 
        
        %Some output variables
        epschi.eps(ii,jj) = epsilon;
        epschi.spec_batch(ii,jj,:) = f11;
        epschi.spec_Y_ratio(ii,jj,kInds) = Y;
        epschi.kmin(ii,jj) = kL;
        epschi.kmax(ii,jj) = kU;
        epschi.kB(ii,jj) = kB;
        epschi.kstar(ii,jj) = k_star;
        epschi.num_k(ii,jj) = length(kInds);
        epschi.X_obs(ii,jj) = X_obs;
        epschi.X_lw(ii,jj) = X_lw;
        epschi.X_hw(ii,jj) = X_hw;
        epschi.MAD(ii,jj) = sum(abs(Y-mean(Y)))/length(Y);

    end
    
    %This is a bit hoaky, but if a spectral calculation failed, the
    %spectrum is returned as complex vector of NaNs. The NaNs are okay, but
    %the complex bit is annoying for later analysis, so ditch it
    epschi.spec_measu = real(epschi.spec_measu);
    
    if stop2draw
        fprintf('Estimate %d/%d of this cast\n',ii,nEstimates);
        disp('Paused. F5 or ''dbcont'' to continue...')
        keyboard
        %Clear all figures before continuing
        figs=findobj(0,'type','figure');
        for mm=cell2mat({figs.Number})
            figure(mm); clf;
        end
    end
    
end %loop through chi and epsilon estimates





















end