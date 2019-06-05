
function Sn = micT_noise(fr)
%USE: Sn = micT_noise(fr)
%    fr - 1D vector of frequencies
%    Sn - 2xN array containing vectors with noise spectra for
%         probes MicT1 (first row) and MicT2 (second row)
%
%return noise spectrum in frequency space for two micT temperature
%probes. These are just empirical functions based on fitting to a "quiet"
%section of the temperature records. You'll probably need to redo this
%everytime you get new probes or have a new mission. If you don't have a
%reasonable noise model you risk getting funky batchelor fits.
%
%the processing parameters used in t02_make_noise_model.m should be the
%same as those used in t03_micT_processing.m, otherwise the vectors won't
%line up. The frequency vector is passed to this function only as a sanity
%check: if it doesn't match the frequency vector from the noise model, you
%screwed up somewhere.
%
%There's also a simple theoretical model commented out below that you can
%implement instead if you prefer. the basic model (m=0) is a white-noise
%model (in log-log space). Changing b shifts the curve up and down,
%changing m tilts the curve. I found that the empirical model works miles
%better when estimating the Batchelor fits.

%B.Scheifele 2017-08

%it would be faster to load once in the parent script and pass the
%noise spectra to micT_eps_chi.m as an input, but I couldn't be bothered to
%change this...

load noise_spectra_geotr2015.mat

if ~isequal(fr,fr_noise)
   error('Something is definitely wrong. The two frequency vectors should be the same'); 
end

Sn = S_noise;

end

%If not using the empirical model, you could use this theoretical model, a
%straight line in log-log space. I started with this (GEOTRACES2015 values
%shown, other previously used values given below) but then moved to the
%empirical model, as developed in t02_make_noise_model.m.
%B.Scheifele, 2017-06-16

% %Probe 1 (blue), conservative estimate
% b = -10.7; m = -0.6;
% Sn1=(10^b)*fr.^m;
% %Probe 2 (red), conservative estimate
% b = -10.45; m = -0.6;
% Sn2=(10^b)*fr.^m;
% %Combined
% Sn = [Sn1; Sn2];



%Previously used values by J. Carpenter
%b = -11.5; m = 0; % Baltic Sea values
%b=-10.6; m=1.15; % Lake Kivu values
%b=-7.5; m=-2.2; % North Sea values
