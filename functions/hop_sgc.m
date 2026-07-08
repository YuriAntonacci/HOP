%% HOP - frequency domain bivariate Granger Causality from ISS model
% modified from Barnett function iss_SGC (it is the unconditional spectral GC applied to restricted model)
% NB: it is different than that produced by OIR=hop_fdTE (this one is the Geweke implementation of Barnett SS GC)

% A,C,K,V: SS model parameters
% i1: indexes of driver
% i2: indexes of target

function [gc,Sout,H] = hop_sgc(A,C,K,V,i1,i2,nfft)

%% preliminaries
omega = ((0:nfft-1)/nfft)*pi; % frequencies in range [0,pi]
z = exp(-1i*omega);         % on unit circle in complex plane

%%% Reduced ISS Model
[CR,KR,VR]=hop_submodel(A,C,K,V,[i1 i2]);
Mr1=length(i1); Mr2=length(i2); Mr=Mr1+Mr2; j1=1:Mr1; j2=Mr1+1:Mr;% indexes of the two blocks inside the reduced process

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(CR); assert(m1 == m);
[m1,n1] = size(KR); assert(n1 == n && m1 == m);
[n1,n2] = size(VR); assert(n1 == n && n2 == n);
z = z(:);

j2 = j2(:)'; assert(all(j2 >=1 & j2 <= n));
j1 = j1(:)'; assert(all(j1 >=1 & j1 <= n));
assert(isempty(intersect(j2,j1)));
assert(length(j1)+length(j2)==n); 
%aggiunta Luca: devono entrare reduced models dove si usano tutti gli indici in analisi unconditioned

i3  = 1:n; i3([j2 j1]) = [];
% i13 = [j2 i3];
i23 = [j1 i3];
% i1R = 1:length(j2);

%% computation
gc = nan(nfft,1);

H = hop_issMA(A,CR,KR,z); %H=iss_MA(A,C,K,z);
VSQRT = chol(VR,'lower');
PVSQRT = chol(hop_parcovar(VR,i23,j2),'lower');

Sout=nan*ones(n,n,nfft);
for k = 1:nfft
    HV   = H(:,:,k)*VSQRT;
    S    = HV*HV'; % CPSD (eq. 6)
    S11  = S(j2,j2);
    HV12 = H(j2,j1,k)*PVSQRT;
    gc(k) = log(real(det(S11))) - log(real(det(S11-HV12*HV12'))); % eq. 11
    Sout(:,:,k)=S;
end

% Note: in theory we shouldn't have to take the real part of the determinants of
% the (Hermitian, positive-definite) matrices in the calculation of the f(k),
% since they should be real. However, Matlab's det() will sometimes return a
% tiny imaginary part.
