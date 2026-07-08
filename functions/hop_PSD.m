%% Power Spectral density from ISS model

% A,C,K,V - parameters of the ISS model to analyze
% Fs - sampling frequency
% nfft - number of points on the frequency axis

% S - Power Spectral Density matrix
% freq - vector of frequencies (nfft points)

function [S,freq] = hop_PSD(A,C,K,V,Fs,nfft)

narginchk(4,6)
if nargin<6, nfft=512; end
if nargin<5, Fs=1; end

M=size(V,1);

freq = (0:nfft-1)*(Fs/(2*nfft)); % frequency axis
s = exp(1j*2*pi*freq/Fs); % vector of complex exponentials
H = hop_issMA(A,C,K,s); % transfer function from ISS parameters (SSGC toolbox Barnett)

S=zeros(M,M,nfft); % Spectral Matrix S(w)
for n=1:nfft % at each frequency
    S(:,:,n)  = H(:,:,n)*V*H(:,:,n)'; %PSD                  
end

end



