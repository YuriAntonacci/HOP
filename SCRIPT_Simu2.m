%% Simulation of PIRD/causal PIRD with N=2 sources
% reproduces Simulation 2 of Faes et. al submitted to IEEE-TBME [2026]
clear; close all; clc;

% addpath([pwd, '/functions/']);

Fs=1; nfft=512;
M=3;
N=M-1;

%% SIMULATION setting
%%% 3 simulations with 3 AR processes: each of the simulation selects a 
% process as target and the other 2 as sources

M=3; Mv=[1 1 1];
par.poles{1}=([0.7 0.1; 0.9 0.35]); % Self-oscillation 1 Oscillations RR
par.poles{2}=([0.7 0.1;]); % Self-oscillation 2
par.poles{3}=([]); % Self-oscillations 3
ncoeff=20;
b1=fir1(ncoeff,2*0.2,'low');
b2=fir1(ncoeff,2*0.2,'high');
par.coup=[2 3 1 1];
c=0.4;
for k=1:ncoeff+1
    par.coup=[par.coup; 1 3 k 1*c*b1(k)];
    par.coup=[par.coup; 1 2 k 1*(1-c)*b2(k)];
end

% par.coup=[par.coup; 1 3 1 0.4];
% par.coup=[par.coup; 1 2 3 0.4];
par.Su=[2 0.5 2]; %variance of innovation processes

[Am,Su]=theoreticalVAR(M,par); % % VAR parameters
Am=Am'; Su=Su';

%%

freq = (0:nfft-1)*(Fs/(2*nfft));

load("Mobius2.mat");
is = nchoosek(1:3, 2); %indexes of sources for the 3 simulations
it = 3:-1:1; %index of target for the 3 simulations

values = cell(3, 10);

for simulation = 1:3

    pird_out = hop_sspird(Am, Su, is(simulation, :), it(simulation), PID, Fs, nfft);
    gc_out = hop_pdgc(Am, Su, is(simulation, :), it(simulation), PID, Fs, nfft);
    
    % The function 'hop_sspird' takes as input the parameters of a VAR model,
    % the innovation covariance matrix, the indices of the source and target
    % variables, the PID lattice structure, and the coarse-graining scheme.
    %
    % It computes the power spectral density (PSD) and, from it, derives the
    % mutual information rate (MIR) and the transfer entropy (TE) between each
    % subset of source variables and the target.
    %
    % For linear VAR processes, TE is equivalent to spectral Granger causality
    % (GC), according to the relation TE = GC / 2.
    %
    % The function then evaluates the remaining redundancy terms of the PID
    % lattice by applying the minimum mutual information criterion between
    % source subsets and the target. The spectral values of the PID atoms are
    % obtained through the Fast Moebius Transform.
    %
    % Finally, the spectral atoms can be integrated to recover the corresponding
    % time-domain measures, or grouped to obtain a coarse-grained representation
    % of the system interactions. This coarse-grained representation can also be
    % mapped back to the time domain through spectral integration.
    %
    % The function 'hop_pdgc' follows the same computational procedure as
    % 'hop_sspird'; however, instead of decomposing the mutual information rate
    % or transfer entropy into PID atoms, it decomposes the spectral Granger
    % causality measures, computed through the function 'hop_sgc'.

    
    values{simulation, 1} = sum(pird_out.cgmir)';
    values{simulation, 2} = [sum(pird_out.cgte)', sum(gc_out.cggc)'/2];
    values{simulation, 3} = pird_out.cgmir(1, :)';
    values{simulation, 4} = pird_out.cgmir(2, :)';
    values{simulation, 5} = pird_out.cgmir(3, :)';
    values{simulation, 6} = pird_out.cgmir(4, :)';
    values{simulation, 7} = [pird_out.cgte(1, :)', gc_out.cggc(1, :)'/2];
    values{simulation, 8} = [pird_out.cgte(2, :)', gc_out.cggc(2, :)'/2];
    values{simulation, 9} = [pird_out.cgte(3, :)', gc_out.cggc(3, :)'/2];
    values{simulation, 10} = [pird_out.cgte(4, :)', gc_out.cggc(4, :)'/2];

end

%% figures

% colNames = {'$\mathbf{MIR}$','$\frac{1}{2}f_{X \rightarrow Y}(f)$','$\mathbf{MIR: Un_1}$',...
%     '$\mathbf{MIR: Un_2}$','$\mathbf{MIR: Syn}$','$\mathbf{MIR: Red}$',...
%     '$\mathbf{GC/2: Un_1}$','$\mathbf{GC/2: Un_2}$','$\mathbf{GC/2: Syn}$',...
%     '$\mathbf{GC/2: Red}$'}';
rowNames = {'$\mathbf{Y = S_3}$','$\mathbf{Y = S_2}$','$\mathbf{Y = S_1}$'};

xLabelText = '$f$ [Hz]';
yLabelText = '[nats]';
font_name = 'Times New Roman';
font_size_big = 18;
font_size_small = 14;
colors = [181 87 181; 24 143 102; 31 160 223; 254 92 92; 0 0 0]/255;
color_order = [5, 5, 1, 2, 3, 4, 1, 2, 3, 4];

figWidth  = 1800;                % pixels
figHeight = round(0.3*figWidth);

figure('Color','w',...
    'Position',[100 100 figWidth figHeight]);

try 

    ha = tight_subplot(3,10,...
        [0.0075 0.0075],...   % gaps [vertical horizontal]
        [0.12 0.08],...   % lower/upper margins
        [0.06 0.02]);     % left/right margins
    
    mm_par = [values{:}]; mm_par = mm_par(:); min_par = min(mm_par); 
    max_par = max(mm_par); 
    
    for r = 1:3
        for c = 1:10
    
            idx = (r-1)*10 + c;
    
            axes(ha(idx));
    
            y = values{r, c};
            plot(freq,y(:, 1),'LineWidth',1.5, 'Color', colors(color_order(c), :));
            ylim([min_par*1.05, max_par*1.05]);
            xlim([-0.01 0.51]);
    
            box on
            grid on
    
            if r == 1
                % title(colNames{c}, ...
                %     'Interpreter','latex', ...
                %     'FontName','', ...
                %     'FontWeight','normal', ...
                %     'FontSize', font_size_big);
            end
    
            if c == 1
                ylabel({rowNames{r}; yLabelText}, ...
                    'FontWeight','bold', 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'YTickLabel',[], 'FontName', font_name, 'FontSize', font_size_small);
            end
    
            if r == 3
                xlabel(xLabelText, 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'XTickLabel',[], 'FontSize', font_size_small);
            end
    
        end
    
    end

catch

    mm_par = [values{:}]; mm_par = mm_par(:); min_par = min(mm_par); 
    max_par = max(mm_par); 

    for r = 1:3
        for c = 1:10

            idx = (r-1)*10 + c;
            subplot(3, 10, idx)

            y = values{r, c};
            plot(freq,y(:, 1),'LineWidth',1.5, 'Color', colors(color_order(c), :));
            ylim([min_par*1.05, max_par*1.05]);
            xlim([-0.01 0.51]);

            box on
            grid on

            % if r == 1
            %     title(colNames{c}, ...
            %         'FontWeight','bold', ...
            %         'Interpreter','latex', 'FontSize', font_size_big);
            % end

            if c == 1
                ylabel({rowNames{r}; yLabelText}, ...
                    'FontWeight','bold', 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'YTickLabel',[], 'FontName', font_name, 'FontSize', font_size_small);
            end

            if r == 3
                xlabel(xLabelText, 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'XTickLabel',[], 'FontSize', font_size_small);
            end

        end

    end

end

%%

figure('Color','w',...
    'Position',[100 100 figWidth figHeight]);

try 

    ha = tight_subplot(3,10,...
        [0.0075 0.0075],...   % gaps [vertical horizontal]
        [0.12 0.08],...   % lower/upper margins
        [0.06 0.02]);     % left/right margins
    
    mm_par = [values{:}]; mm_par = mm_par(:); min_par = min(mm_par); 
    max_par = max(mm_par); 
    
    for r = 1:3
        for c = 1:10
    
            idx = (r-1)*10 + c;
    
            axes(ha(idx));
    
            y = values{r, c};
            try
                plot(freq,y(:, 2),'LineWidth',1.5, 'Color', colors(color_order(c), :));
            catch
                plot(freq,y,'LineWidth',1.5, 'Color', colors(color_order(c), :));
            end
            ylim([min_par*1.05, max_par*1.05]);
            xlim([-0.01 0.51]);
    
            box on
            grid on
    
            % if r == 1
            %     title(colNames{c}, ...
            %         'FontWeight','bold', ...
            %         'Interpreter','latex', 'FontSize', font_size_big);
            % end
    
            if c == 1
                ylabel({rowNames{r}; yLabelText}, ...
                    'FontWeight','bold', 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'YTickLabel',[], 'FontName', font_name, 'FontSize', font_size_small);
            end
    
            if r == 3
                xlabel(xLabelText, 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'XTickLabel',[], 'FontSize', font_size_small);
            end
    
        end
    
    end

catch

    mm_par = [values{:}]; mm_par = mm_par(:); min_par = min(mm_par); 
    max_par = max(mm_par); 

    for r = 1:3
        for c = 1:10

            idx = (r-1)*10 + c;
            subplot(3, 10, idx)

            y = values{r, c};
            try
                plot(freq,y(:, 2),'LineWidth',1.5, 'Color', colors(color_order(c), :));
            catch
                plot(freq,y,'LineWidth',1.5, 'Color', colors(color_order(c), :));
            end
            ylim([min_par*1.05, max_par*1.05]);
            xlim([-0.01 0.51]);

            box on
            grid on

            % if r == 1
            %     title(colNames{c}, ...
            %         'FontWeight','bold', ...
            %         'Interpreter','latex', 'FontSize', font_size_big);
            % end

            if c == 1
                ylabel({rowNames{r}; yLabelText}, ...
                    'FontWeight','bold', 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'YTickLabel',[], 'FontName', font_name, 'FontSize', font_size_small);
            end

            if r == 3
                xlabel(xLabelText, 'Interpreter','latex', 'FontSize', font_size_big);
            else
                set(gca,'XTickLabel',[], 'FontSize', font_size_small);
            end

        end

    end

end