%% frequency domain bivariate Transfer Entropy from reduced ISS model
% identical (except for 0.5 factor) to oir_fdGC of the OIR toolbox ("MIR" approach of Faes IEEE-T-SigProc)

% A,CR,KR,VR - parameters of the ISS reduced model to analyze
% j1,j2 - vector indexes of the two blocks for which to compute the bivariate measures
% j1 and j2 must cover the whole dimension of the reduced model: length(j1)+length(j2)=Mr, [j1 j2]=1:Mr

function [f12,f1_2,f2_1,f1o2,S,freq] = hop_fdTE(A,CR,KR,VR,j1,j2,Fs,nfft)

narginchk(6,8)
if nargin<8, nfft=512; end
if nargin<7, Fs=1; end

Mr=length(j1)+length(j2);
assert(Mr==size(CR,1));
% assert(sum([j1 j2]==1:Mr)==Mr); %commentata per il toolbox SGC..
assert( sum([j1 j2]) == sum(1:Mr) ); %nuova per il toolbox SGC (dovrebbe essere equivalente)

bVR{1,1}=VR(j1,j1); bVR{1,2}=VR(j1,j2);
bVR{2,1}=VR(j2,j1); bVR{2,2}=VR(j2,j2);

freq = (0:nfft-1)*(Fs/(2*nfft)); % frequency axis
s = exp(1j*2*pi*freq/Fs); % vector of complex exponentials
H = hop_issMA(A,CR,KR,s); % transfer function from ISS parameters (SSGC toolbox Barnett)

S=zeros(Mr,Mr,nfft); % Spectral Matrix S(w)
bH=cell(2,2,nfft); 
bS=cell(2,2,nfft); 
f12=nan*ones(nfft,1);
f1_2=nan*ones(nfft,1); f2_1=f1_2; % Geweke spectral causality from sources to destinations
for n=1:nfft % at each frequency
    S(:,:,n)  = H(:,:,n)*VR*H(:,:,n)'; %PSD
    
    %%% extraction of blocks indexes
    bH{1,1,n}=H(j1,j1,n); bH{1,2,n}=H(j1,j2,n);
    bH{2,1,n}=H(j2,j1,n); bH{2,2,n}=H(j2,j2,n);  
    bS{1,1,n}=S(j1,j1,n); bS{1,2,n}=S(j1,j2,n);
    bS{2,1,n}=S(j2,j1,n); bS{2,2,n}=S(j2,j2,n);
    
    %instterm=bH{2,2,n}*bVR{2,1}*bH{2,1,n}' + bH{2,1,n}*bVR{1,2}*bH{2,2,n}';
    %f(n) = 0.5*log( abs(det(bS{2,2,n})) / abs(det(bS{2,2,n}-bH{2,1,n}*bVR{1,1}*bH{2,1,n}'- instterm)) );
    f2_1(n) = 0.5*log( abs(det(bS{1,1,n})) / abs(det(bH{1,1,n}*bVR{1,1}*bH{1,1,n}')));
    f1_2(n) = 0.5*log( abs(det(bS{2,2,n})) / abs(det(bH{2,2,n}*bVR{2,2}*bH{2,2,n}')));
    
    f12(n)=0.5*log( abs(det(bS{1,1,n}))*abs(det(bS{2,2,n})) / abs(det(S([j1 j2],[j1 j2],n))) );
                   
end
f1o2=f12-f1_2-f2_1;

end



