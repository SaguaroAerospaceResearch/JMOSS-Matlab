% JMOSS V4 Demo File
% This script demonstrates how to use the basic features of the JMOSS Air
% Data Calibration algorithm.
%
% Written by Juan Jurado, Air Force Institute of Technology, 2018

clear; close all; clc;

%% Load sample data and feed to algorithm
% In this example, we will load two test points simultaneously to exercise
% the neural net training features with more than a single test point.
load sampleData.mat;
%#ok<*SAGROW>
for ii = 1:length(data)
    TestPoints(ii).Pt = data(ii).Pt; % [Psi]
    TestPoints(ii).Ps = data(ii).Ps; % [Psi]
    TestPoints(ii).Tic = data(ii).Tic; % [K]
    TestPoints(ii).alphai = data(ii).alphai; % [rad]
    TestPoints(ii).betai = data(ii).betai; % [rad]
    TestPoints(ii).Vn = data(ii).Vn; % [ft/s]
    TestPoints(ii).Ve= data(ii).Ve;% [ft/s]
    TestPoints(ii).Vd = data(ii).Vd; % [ft/s]
    TestPoints(ii).hg = data(ii).hg; % [ft MSL]
    TestPoints(ii).roll = data(ii).roll; % [rad]
    TestPoints(ii).pitch = data(ii).pitch; % [rad]
    TestPoints(ii).yaw = data(ii).yaw; % [rad]
    TestPoints(ii).time = data(ii).time; % [seconds]
end

%% Feed TestPoints to JMOSSV4
[MODEL,EKF] = JMOSSV4(TestPoints);
% Other example usage:
% [MODEL,EKF] = JMOSSV4(TestPoints(1)); % Only feed one of the TestPoints at a time
% [MODEL,EKF] = JMOSSV4(TestPoints,'alpha',0.1); % Specify a confidence level of 90% (default 95%)
% [MODEL,EKF] = JMOSSV4(TestPoints,'rollLimit',15); % Specify a AOB limit for excluding deltaP results (default is AOB>10 deg)

%% Process results for plotting
% Raw EKF results
ekfmach = cell2mat({EKF.mach}'); % Raw Extended Kalman Filter (EKF) mach vector (all test points)
ekfdPp_Ps = cell2mat({EKF.dPp_Ps}');% Raw EKF dPp_Ps vector (all test points)
turnPts = cell2mat({EKF.deweightedPts}');  % Index of excluded points due to AOB

% GRP Neural Net (Smoothed Output) 
grpMach = MODEL.mach; % Smoothed GRP mach vector (all test points)
grpdPp_Ps = MODEL.dPp_Ps; % Smoothed GRP dPp_Ps curve (all test points)
pb = MODEL.PredictionBand; % 95% Prediction Band

%% Plot main results
% Display additional data from TestPoint 1.
TP = 1;
fontSize = 16;
WnHat = mean(EKF(TP).bkdPass(:,2)); % North wind [fps]
WeHat = mean(EKF(TP).bkdPass(:,3)); % East wind [fps]
WdHat = mean(EKF(TP).bkdPass(:,4)); % Down wind [fps]
rRange = EKF(TP).bkdPass([end 1],5); % Range of Temperature recovery factor values [unitless]

stateStrs = {sprintf('\\bf{Additional Parameters (Test Point %0.0f)}',TP)
             sprintf('$\\hat{V}_{W_n}$ = %0.5f ft/s',WnHat);
             sprintf('$\\hat{V}_{W_e}$ = %0.5f ft/s',WeHat);
             sprintf('$\\hat{V}_{W_d}$ = %0.5f ft/s',WdHat);
             sprintf('$\\hat{r} \\in [%0.3f,%0.3f]$',rRange)};
             
figure;
h(1) = plot(ekfmach(~turnPts),ekfdPp_Ps(~turnPts),'bo','MarkerFaceColor','b',...
    'MarkerSize',3); hold on;
h(2) = plot(ekfmach(turnPts),ekfdPp_Ps(turnPts),'go','MarkerFaceColor','g',...
    'MarkerSize',3);
h(3) = plot(grpMach,grpdPp_Ps,'r-','LineWidth',3);
h(4:5) = plot(grpMach,pb,'r--','LineWidth',2);
set(gca,'FontSize',fontSize,'FontName','helvetica');
t = text(0,0,stateStrs,'Interpreter','latex','EdgeColor','k',...
    'BackgroundColor',[1 1 1],'FontSize',fontSize,'Units','normalized');
set(t,'Position',[0.01,0.77,0])

xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('SPE, $\Delta P_p/P_s$','Interpreter','latex',...
    'FontSize',fontSize);
axis tight;
l = legend(h([1 2 3 4]),'Raw EKF Ouput','Excluded Turn Data','GRP','$95$\% Pred. Band',...
    'Location','NorthWest');
set(l,'Interpreter','latex','FontSize',fontSize);
grid minor;
set(gcf,'Position',[0 0 1.5*800 800]);
title('\textbf{Static Position Error, JMOSS V4 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)

%% AoA Effects
% We can look at each test point individually to see how they might differ
% from run to run. Be on the lookout for dPp_Ps vs. machCI changes due to 
% chanes in w/delta, which can be diagnosed by plotting against AoA.  In
% this plot, we'll go ahead and not plot the turn data so it doesn't clutter
% the plot. We can animate the 3D plot so get a good look at how AoA affects
% dPp_Ps vs. MachIC.
animate = true;
figure; hold on;
nPoints = length(data);
colors = jet(nPoints);
names = cell(nPoints,1);
for ii = 1:nPoints
    subSamp = EKF(ii).subSample; % Get the indices of decimated points 
    c = colors(ii,:);
    mach = EKF(ii).mach;
    dPpPs = EKF(ii).dPp_Ps;
    alpha = TestPoints(ii).alphai(subSamp);
    turn = EKF(ii).deweightedPts;
    h1(ii) = plot3(mach(~turn),dPpPs(~turn),alpha(~turn),'o','Color',c,'MarkerFaceColor',c,'MarkerSize',3);
    names{ii} = sprintf('TestPoint %0.0f',ii);
end
xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('SPE, $\Delta P_p/P_s$','Interpreter','latex',...
    'FontSize',fontSize);
zlabel('AoA [rad]');
axis tight;
legend(h1,names)
grid minor;
set(gcf,'Position',[0 0 1.5*800 800]);
title('\textbf{Static Position Error: AoA Effects, JMOSS V4 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)
view(3);
if animate
    M = 200;
    az = linspace(40,-20,M);
    el = 30;
    for ii = 1:M
        view(az(ii),el);
        drawnow;
    end
end
%% State Diagnostics
% Here we can look at the EKF state history on the first and second (final)
% pass. This can be used to look at EKF estimates of 3D wind, total temp
% recovery factor (r) and Pa linearization point (P0). We will focus on the
% first TestPoint just to see an example.
TP = 1;
stateTitles = {'$\Delta P_p/P_s$','$V_{W_N}$','$V_{W_E}$','$V_{W_D}$','$r$','$P_0$'};
N = length(stateTitles);
figure;
for ii = N:-1:1
    subplot(N,1,ii);
    plot(EKF(TP).mach,EKF(TP).fwdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 6
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end
titleStr = sprintf('\\bf{JMOSS EKF: Forward Pass (Test Point %0.0f)}',TP);
title(titleStr,'Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 800 1600]);

figure;
for ii = N:-1:1
    subplot(N,1,ii);
    plot(EKF(TP).mach,EKF(TP).bkdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 6
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end
titleStr = sprintf('\\bf{JMOSS EKF: Backward Pass (Test Point %0.0f)}',TP);
title(titleStr,'Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 800 1600]);








