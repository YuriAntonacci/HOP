%% Computation of MIR and of its time-domain decomposition and frequency-domain expansion for two blocks of processes
%  identical (except to 0.5 factor) to oir_mir of the OIR toolbox

function out=hop_mir(A,C,K,Su,Mv,i_1,i_2,Fs,nfft)

narginchk(7,10)
if nargin<9, nfft=512; end
if nargin<8, Fs=1; end %default non compute time domain measures through submodels

% indexes of the two blocks to analyze inside the Q time series
[i1,i2]=hop_subindexes(Mv,i_1,i_2);

% reduced model with the two blocks to analyze - ired=[i1 i2]
[CR,KR,VR,rep]=hop_submodel(A,C,K,Su,[i1 i2]); % Eqs. (27-28)
if rep<0
    out.rep=-1;
    return
end
out.VR=VR;

% indexes of the two blocks inside the reduced process
Mr1=length(i1); Mr2=length(i2); Mr=Mr1+Mr2;
j1=1:Mr1; j2=Mr1+1:Mr;

% FREQUENCY DOMAIN TE
[t12,t1_2,t2_1,t1o2,S,freq] = hop_fdTE(A,CR,KR,VR,j1,j2,Fs,nfft);

T2_1=mean(t2_1); %TE from i2 to i1
T1_2=mean(t1_2); %TE from i1 to i2
T1o2=mean(t1o2); 
I12=mean(t12);

out.I12=I12;
out.T1_2=T1_2;
out.T2_1=T2_1;
out.I1o2=T1o2;
out.t12=t12;
out.t1_2=t1_2;
out.t2_1=t2_1;
out.t1o2=t1o2;
out.S=S;
out.freq=freq;
out.rep=1;


end