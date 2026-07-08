%% Application of PDGC with N=2 sources to a representative set of cardiovascular time series
clear; close all; clc;

nfft=512;
LF=[0.03 0.15];
HF=[0.15 0.4];

%% open data, preprocess
alldata=readtable('series.csv'); alldata=alldata{:,:};
rr=alldata(:,1);
sap=alldata(:,2);
map=alldata(:,3);
cbfv=alldata(:,4);
resp=alldata(:,5);

data=[resp sap rr];
[Ns,M]=size(data);
N=M-1;

Fs=1/mean(rr); %sampling frequency based on mean heart period


%%% normalization of series
X=nan*ones(Ns,M);
for m=1:M 
    X(:,m)=(data(:,m)-mean(data(:,m)))./std(data(:,m)); %Set 0 mean and 1 stdev
end

%%% arrange bands
nLF=round(LF.*nfft/(Fs/2));
nHF=round(HF.*nfft/(Fs/2));

%% analysis
load('Mobius2.mat');
L=size(PID.ac,2);

%%% VAR  model order selection and identification
[pottaic,pottbic] = hop_VARmos(X',12,0); %order selection
p=pottbic;

[Am,Su]=hop_VARid(X',p,0); %identification

%%% High-order GC analysis
is=[1 2]; %indexes of sources
it=3; %index of target
out=hop_pdgc(Am,Su,is,it,PID,Fs,nfft);

cgGC_LF=sum(out.cggc(:,nLF(1):nLF(2)-1),2)./nfft;
cgGC_HF=sum(out.cggc(:,nHF(1):nHF(2)-1),2)./nfft;
aGC_LF=sum(out.agc(:,nLF(1):nLF(2)-1),2)./nfft;
aGC_HF=sum(out.agc(:,nHF(1):nHF(2)-1),2)./nfft; 
rGC_LF=sum(out.rgc(:,nLF(1):nLF(2)-1),2)./nfft;
rGC_HF=sum(out.rgc(:,nHF(1):nHF(2)-1),2)./nfft;

      

%% plots
lw=1.5;

Rlab{1}='{X_{1}\rightarrow Y}';
Rlab{2}='{X_{2}\rightarrow Y}';
Rlab{3}='{X_{\{1,2\}}\rightarrow Y}';
Rlab{4}='{X_{\{1\}\{2\}}\rightarrow Y}';

ymin=min([0 min(out.rgc(:))]); ymax=max(out.rgc(:));
figure('name','Redundant GC')
for i=1:N+2
    subplot(2,2,i);
    plot(out.freq,out.rgc(i,:),'color',[0.4 0.4 0.4],'linewidth',lw); %hold on
    xlim([0 max(out.freq)]); ylim([ymin ymax]);
    xlabel('f')
    ylabel(['f_' Rlab{i}])
    title(['F_' Rlab{i} '=' num2str(out.rGC(i))])
end

col=[182 88 182; 25 143 103; 32 160 224; 255 92 92]./255;
ymin=min([0 min(out.cggc(:))]); ymax=max(out.cggc(:));
figure('name','PDGC')
for i=1:N+2
    subplot(2,2,i);
    plot(out.freq,out.cggc(i,:),'color',col(i,:),'linewidth',lw); %hold on
    xlim([0 max(out.freq)]); ylim([ymin ymax]);
    xlabel('f')
    ylabel([out.cglab(i,:) '(f)'])
    title([out.cglab(i,:) '=' num2str(out.cgGC(i))])
end

