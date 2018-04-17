% Example use of JMOSS V3.1 Pitot-static calibration algorithm
%
% Authors: Juan Jurado and Clark McGehee, Copyright 2018
clc; clear; close all;

%% Load Data and Construct Test Point Structure
% Note: You can load multiple test points into the JMOSS algorithm by
% adding to the TestPoint structure (i.e. TestPoint(1).Pt, TestPoint(2).Pt,
% etc...). The final smoothing will occur on the combined results from all
% test points loaded in this manner.

data = load('SampleData.mat'); % UNITS:
TestPoints(1).Pt = data.Pt; % [Psi]
TestPoints(1).Ps = data.Ps; % [Psi]
TestPoints(1).Tic = data.Tic; % [K]
TestPoints(1).alphai = data.alphai; % [rad]
TestPoints(1).betai = data.betai; % [rad]
TestPoints(1).Vn = data.Vn; % [ft/s]
TestPoints(1).Ve= data.Ve;% [ft/s]
TestPoints(1).Vd = data.Vd; % [ft/s]
TestPoints(1).hg = data.hg; % [ft MSL]
TestPoints(1).roll = data.roll; % [rad]
TestPoints(1).pitch = data.pitch; % [rad]
TestPoints(1).yaw = data.yaw; % [rad]
TestPoints(1).time = data.time; % [seconds]
%% Feed TestPoint structure to JMOSSV31
[ASM,EKFOutput] = JMOSSV31(TestPoints); % Uses all default settings
% Other example usage:
% [ASM,EKFOutput] = JMOSSV31(TestPoints,'alpha',0.1);
% [ASM,EKFOutput] = JMOSSV31(TestPoints,'freezeStates',false);
% [ASM,EKFOutput] = JMOSSV31(TestPoints,'knots',.9:.01:1.0,'alpha',0.01,'freezeStates',false);
% [ASM,EKFOutput] = JMOSSV31(TestPoints,'maxKnots',10);

%% Process results for plotting
ekfmach = cell2mat({EKFOutput.mach}'); % Raw EKF mach vector (all test points)
ekfdPp_Ps = cell2mat({EKFOutput.dPp_Ps}');% Raw EKF dPp_Ps vector (all test points)
asmmach = ASM.mach; % Smoothed ASM mach vector (all test points)
asmdPp_Ps = ASM.dPp_Ps; % Smoothed ASM dPp_Ps curve (all test points)
pb = ASM.PredictionBand; % 95% Prediciton Band (i.e. Tolerance Interval) 

%% Plot Results Including Wind and Kt estiamtes
fontSize = 16;
WnHat = mean(EKFOutput(1).bkdPass(:,2));
WeHat = mean(EKFOutput(1).bkdPass(:,3));
WdHat = mean(EKFOutput(1).bkdPass(:,4));
KtHat = mean(EKFOutput(1).bkdPass(:,5));

stateStrs = {'\bf{Additional Parameters:}';
             sprintf('$\\hat{V}_{W_n}$ = %0.5f ft/s',WnHat);
             sprintf('$\\hat{V}_{W_e}$ = %0.5f ft/s',WeHat);
             sprintf('$\\hat{V}_{W_d}$ = %0.5f ft/s',WdHat);
             sprintf('$\\hat{K}_t$ = %0.5f',KtHat)};

figure;
h(1) = plot(ekfmach,ekfdPp_Ps,'bo','MarkerFaceColor','b',...
    'MarkerSize',3); hold on;
h(2) = plot(asmmach,asmdPp_Ps,'r-','LineWidth',3);
h(3:4) = plot(asmmach,pb,'r--','LineWidth',2);

set(gca,'FontSize',fontSize,'FontName','helvetica');
text(0.55,0.04,stateStrs,'Interpreter','latex','EdgeColor','k',...
    'BackgroundColor',[1 1 1],'FontSize',fontSize);

xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('SPE, $\Delta P_p/P_s$','Interpreter','latex',...
    'FontSize',fontSize);
axis tight;
l = legend(h([1 2 3]),'Raw EKF Ouput','ASM','$95$\% Pred. Band',...
    'Location','NorthWest');
set(l,'Interpreter','latex','FontSize',fontSize);
grid minor;
set(gcf,'Position',[0 0 1.5*800 800]);
title('\textbf{Static Position Error, JMOSS Algorithm Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)

%% Diagnostics
% Here we can look at the EKF state history on the first and second (final)
% pass. This can be used to look at EKF estimates of 3D wind, Kt, and
% non-standard pressure correction.
stateTitles = {'$\Delta P_p/P_s$','$V_{W_N}$','$V_{W_E}$',...
    '$V_{W_D}$','$K_t$','$\delta P_0$'};
figure;
for ii = 6:-1:1
    subplot(6,1,ii);
    plot(EKFOutput(1).mach,EKFOutput(1).fwdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 6
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end

title('\bf{EKF: First Pass}','Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 800 1600]);

figure;
for ii = 6:-1:1
    subplot(6,1,ii);
    plot(EKFOutput(1).mach,EKFOutput(1).bkdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 6
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end

title('\bf{EKF: Second Pass}','Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 800 1600]);






