%% Simulation of PDGC with N=3 sources
% reproduces Simulation 1 of Faes et. al submitted to IEEE-TBME [2026]
clear; close all; clc;

Fs=1; nfft=512;
M=4;
N=M-1;
is=[1 2 3]; %indexes of sources
it=4; %index of target

%% SIMULATION setting
%%% simulation with 2 sources and 1 target
par.poles{4}=([]); % Self-oscillation 4 (absent)
par.poles{2}=([0.9 0.3]); % Self-oscillation 2
par.poles{3}=([0.9 0.4]); % Self-oscillations 3
par.poles{1}=([0.94 0.1]); % Self-oscillation 1
par.coup=[2 4 1 1; 3 4 1 1; 1 2 1 1];
par.Su=[2 1.5 0.75 1]; %variance of innovation processes

[Am,Su,Ak,z]=theoreticalVAR(M,par); %generates VAR parameters from poles and coupling coeffs
Am=Am'; Su=Su';


%% analysis
%%%%% power spectral densities
[A,C,K,~] = hop_ar2iss(Am);
PSD = hop_PSD(A,C,K,Su,Fs,nfft); %PSD matrix
P=nan*ones(nfft,M); %individual PSDs
for m=1:M
    P(:,m)=abs(squeeze(PSD(m,m,:)));
end

%%%%% PDGC analysis
load('Mobius3.mat') %yields the struct PID

% GC functions (obtained thorugh ISS-based tool, Barnett)
out=hop_pdgc(Am,Su,is,it,PID,Fs,nfft);
rgc=out.rgc; rGC=out.rGC; %redundant GC (frequency and time domain)
agc=out.agc; aGC=out.aGC; %atomic GC (frequency and time domain)
cggc=out.cggc; cgGC=out.cgGC; %coarse-grained GC (frequency and time domain)
cglab=out.cglab; %coarse-graining labels
freq=out.freq; %frequency axis


%% plots
ac=PID.ac'; %antichains
L=length(ac); %number of antichains

lw=1.5;
col_rgc=colormap(jet(L)); %col_rgc=colormap(winter(L));
col_agc=colormap(hsv(L)); %col_agc=colormap(cool(L));
close

%%% prepare strings for labels/titles
convfact{1}=1; convfact{2}=[10 1]'; convfact{3}=[100 10 1]'; convfact{4}=[1000 100 10 1]'; %(max 4 sources)
acstr=cell(L,1);
for i=1:L
    if ~iscell(ac{i})
        acstr{i}=int2str(is(ac{i})*convfact{length(ac{i})});
    else
        newstr='';
        for j=1:length(ac{i})
            newstr=strcat(newstr,'-',int2str(is(ac{i}{j})*convfact{length(ac{i}{j})}));  
        end
        acstr{i}=newstr(2:end);
    end
end


%%% Figs
h1=figure('name','Power spectral densities');
h2=figure('name','Redundant spectral GCs');
h3=figure('name','Atomic spectral GCs');
h4=figure('name','Coarse-grained spectral GCs');


figure(h1)
for m=1:M
    subplot(2,2,m);
    plot(freq,P(:,m),'color',[0.5 0.5 0.5],'linewidth',lw); %hold on;
    xlabel('f')
    stringa=strcat('P_{',num2str(m),'}(f)');
    title(stringa)
end

Rmin=min([0 min(rgc(:))]); Rmax=max(rgc(:));
Amin=min([0 min(agc(:))]); Amax=max(agc(:));
for il=1:L
    figure(h2);subplot(3,6,il);
    plot(freq,rgc(il,:),'color',col_rgc(il,:),'linewidth',lw); %hold on;
    xlim([0 max(freq)]); ylim([Rmin Rmax]),
    stringa=strcat('R_{',num2str(il),'} (f)'); ylabel(stringa); 
    xlabel('f')
    stringa=strcat('R_{',string(acstr{il}),'\rightarrow',num2str(it),'}=',num2str(0.5*rGC(il)));
    title(stringa)
    
    figure(h3);subplot(3,6,il);
    plot(freq,agc(il,:),'color',col_agc(il,:),'linewidth',lw); %hold on;
    xlim([0 max(freq)]); ylim([Amin Amax]),
    stringa=strcat('A_{',num2str(il),'} (f)'); ylabel(stringa); 
    xlabel('f')
    stringa=strcat('A_{',string(acstr{il}),'\rightarrow',num2str(it),'}=',num2str(0.5*aGC(il))); title(stringa)
end

col_agcN2=colormap(hsv(N+2));
Cmin=min([0 min(cggc(:))]); Cmax=max(cggc(:));
figure(h4)
for i=1:N+2
    subplot(2,3,i);
    plot(freq,cggc(i,:),'color',col_agcN2(i,:),'linewidth',lw); %hold on
    xlim([0 max(freq)]); ylim([Cmin Cmax]);
    xlabel('f')
    ylabel([cglab(i,:) '(f)'])
    title([cglab(i,:) '=' num2str(cgGC(i)) ',' num2str(0.5*cgGC(i))])
end





