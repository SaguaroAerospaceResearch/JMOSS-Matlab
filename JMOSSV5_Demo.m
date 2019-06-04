% JMOSS V5 Demo File
% This script demonstrates how to use the basic features of the JMOSS Air
% Data Calibration algorithm.
%
% Written by Juan Jurado, Air Force Institute of Technology, 2019

clear; close all; clc;

%% Load sample data and feed to algorithm
% In this example, we will load two test points simultaneously to exercise
% the Gaussian Regression Process training features with more than a single test point.
load sampleData.mat;
nPoints = length(data);
TestPoints(nPoints).Pt = [];

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
    
    % NOTE: If you have your own OAT, simply feed it by adding a Ta field
    % into TestPoints. Be sure to NOT include the Tic entry, otherwise
    % JMOSS will ignore Ta and still compute Ta based on Tic.
    % Example: TestPoints(ii).Ta = data(ii).Ta;
end

%% Feed TestPoints to JMOSSV5
[MODEL,EKF] = JMOSSV5(TestPoints);
% Other example usage:
% [MODEL,EKF] = JMOSSV5(TestPoints(1)); % Only feed one of the TestPoints at a time
% [MODEL,EKF] = JMOSSV5(TestPoints,'alpha',0.1); % Specify a confidence level of 90% (default 95%)
% [MODEL,EKF] = JMOSSV5(TestPoints,'rollLimit',15); % Specify a AOB limit for excluding deltaP results (default is AOB>10 deg)

%% Process results for plotting
% Raw EKF results
ekfmach = cell2mat({EKF.mach}'); % Raw Extended Kalman Filter (EKF) mach vector (all test points)
ekfdPp_Ps = cell2mat({EKF.dPp_Ps}');% Raw EKF dPp_Ps vector (all test points)
turnPts = cell2mat({EKF.deweightedPts}');  % Index of excluded points due to AOB

% GRP Regression (Smoothed Output) 
grpMach = MODEL.mach; % Smoothed GRP mach vector (all test points)

grpdPp_Ps = MODEL.dPp_Ps; % Smoothed GRP dPp_Ps curve (all test points)
pb = MODEL.PredictionBand; % 95% prediction band for deltaPp/Ps

grpDeltaH = MODEL.Curves.deltaH; % Smoothed GRP deltaH curve (all test points) [assumes Sea Level]
grpDeltaHPB = MODEL.Curves.deltaHPredBand; % Prediction band for deltaH

grpDeltaV = MODEL.Curves.deltaV; % Smoothed GRP deltaH curve (all test points) [assumes Sea Level]
grpDeltaVPB = MODEL.Curves.deltaVPredBand; % Prediction band for deltaV

grpDeltaM = MODEL.Curves.deltaM; % Smoothed GRP deltaH curve (all test points) [assumes Sea Level]
grpDeltaMPB = MODEL.Curves.deltaMPredBand; % Prediction band for deltaM


%% Plot main results
% Display additional data from TestPoint 1.
TP = 1;
fontSize = 12;
windowSize = 600;

WnHat = mean(EKF(TP).bkdPass(:,2)); % North wind [fps]
WeHat = mean(EKF(TP).bkdPass(:,3)); % East wind [fps]
WdHat = mean(EKF(TP).bkdPass(:,4)); % Down wind [fps]

stateStrs = {sprintf('\\bf{Wind Estimate (Test Point %0.0f)}',TP)
             sprintf('$\\hat{V}_{W_n}$ = %0.5f ft/s',WnHat);
             sprintf('$\\hat{V}_{W_e}$ = %0.5f ft/s',WeHat);
             sprintf('$\\hat{V}_{W_d}$ = %0.5f ft/s',WdHat)};
             
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
set(gcf,'Position',[0 0 1.6*windowSize windowSize]);
title('\textbf{Static Position Error, JMOSS V5 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)

%% Altitude, Airspeed, and Mach number corrections
% JMOSSV5 now computes altitude (deltaH), airspeed (deltaV), and Mach
% number (deltaM) corrections and their associated prediction bands. Let's 
% plot them.
figure; 
subplot(3,1,1);
plot(grpMach,grpDeltaH,'r-','LineWidth',2); hold on;
plot(grpMach,grpDeltaHPB,'r:','LineWidth',2);
grid minor;
title('\textbf{Altitude, Airspeed, and Mach Number Corrections, JMOSS V5 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)
xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('$\Delta H\ [H_c = H_{ic}+\Delta H]$, ft','FontSize',fontSize,'Interpreter','latex');

subplot(3,1,2);
plot(grpMach,grpDeltaV,'r-','LineWidth',2); hold on;
plot(grpMach,grpDeltaVPB,'r:','LineWidth',2);
grid minor;
xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('$\Delta V\ [V_c = V_{ic}+\Delta V]$, kts','FontSize',fontSize,'Interpreter','latex');

subplot(3,1,3);
plot(grpMach,grpDeltaM,'r-','LineWidth',2); hold on;
plot(grpMach,grpDeltaMPB,'r:','LineWidth',2);
grid minor;
xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex'...
    ,'FontSize',fontSize);
ylabel('$\Delta M\ [M_c = M_{ic}+\Delta M]$, M','FontSize',fontSize,'Interpreter','latex');
set(gcf,'Position',[0 0 windowSize 1.6*windowSize])

%% AoA Effects
% We can look at each test point individually to see how they might differ
% from run to run. Be on the lookout for dPp_Ps vs. machCI changes due to 
% chanes in w/delta, which can be diagnosed by plotting against AoA.  In
% this plot, we'll go ahead and not plot the turn data so it doesn't clutter
% the plot.

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
set(gcf,'Position',[0 0 1.6*windowSize windowSize]);
title('\textbf{Static Position Error: AoA Effects, JMOSS V5 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)
view(3);
%% EKF Diagnostics
% Here we can look at the EKF state history on the first and second (final)
% pass. This can be used to look at EKF estimates of 3D wind, total temp
% recovery factor (r) and Pa linearization point (P0). We will focus on the
% first TestPoint just to see an example.
TP = 1;
stateTitles = {'$\Delta P_p/P_s$','$V_{W_N}$','$V_{W_E}$','$V_{W_D}$','$P_0$'};
figure;
for ii = 5:-1:1
    subplot(5,1,ii);
    plot(EKF(TP).mach,EKF(TP).fwdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 5
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end
titleStr = sprintf('\\bf{JMOSS EKF: Forward Pass (Test Point %0.0f)}',TP);
title(titleStr,'Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 windowSize 1.6*windowSize]);

figure;
for ii = 5:-1:1
    subplot(5,1,ii);
    plot(EKF(TP).mach,EKF(TP).bkdPass(:,ii),'b-','LineWidth',2); hold on;
    set(gca,'FontSize',fontSize,'FontName','helvetica');
    grid minor;
    axis tight;
    ylabel(stateTitles{ii},'Interpreter','latex','FontSize',fontSize);
    if ii == 5
        xlabel('Instrument Corrected Mach, $M_{ic}$',...
            'Interpreter','latex','FontSize',fontSize);
    end
end
titleStr = sprintf('\\bf{JMOSS EKF: Backward Pass (Test Point %0.0f)}',TP);
title(titleStr,'Interpreter','latex','FontSize',fontSize+2);
set(gcf,'Position',[0 0 windowSize 1.6*windowSize]);

%% Temperature Diagnostics
% One big factor in obtaining accurate JMOSS results is the accuracy of the
% ambient temperature (or OAT) measurements. By default, JMOSS expects Total
% Temperature (Tic) measurements and will perform a non-linear optimization
% routine to estimate OAT based on Tic, GPS altitude, and Mach number. If
% you have your own OAT measurements, you can substitue them in by
% providing a TestPoint.Ta (instead of Tic) input into the TestPoint
% object. Experimentally, we have noticed that an error of 1 Kelvin in the
% accuracy of OAT (wether it is user-provided or derived by JMOSS) will
% bias the JMOSS results by up to 1.5e-3 dPp/Ps units or approximately 40
% feet, 2 kts, or 0.005 Mach (assuming sea level and M=0.2). This is why it
% is important to supplement JMOSS test points with at least one Tower Fly
% By per w/delta band for your aircraft, in order to "anchor" the resulting
% curves and cancel out any biases from innacurate OAT.

% Let's look at the JMOSS optimizer results for OAT estimation based on the
% TestPoint.Tic data. 

TP = 1; % Let's look at TestPoint 1
if isfield(TestPoints(TP),'Tic')
% If user provided Tic, then we can look at JMOSS OAT estimation
% performance. Otherwise, if the user provided Ta, JMOSS uses it directly
% and there is no OAT estimation step.
    figure; 
    subplot(3,1,1);
    % There are two estimates of OAT we can look at. This first one is the
    % one used by the EKF as part of the deltaPp/Ps estimation process.
    Mic = EKF(TP).mach;
    estimatedOAT1 = EKF.Ta;
    plot(Mic,estimatedOAT1,'b-','LineWidth',2); hold on;
    
    % This second one is the best approximation of OAT based soley on the
    % JMOSS-computed temperature correction factor, Kt. It is by design a
    % bit noisier, but can be used as the temperature calibration for the
    % total temperature sensor.
    Tic = TestPoints(TP).Tic(EKF(TP).subSample);
    Kt = EKF(TP).Kt;
    estimatedOAT2 = Tic./(1+0.2.*Kt.*Mic.^2);
    plot(Mic,estimatedOAT2,'g-','LineWidth',2);
    grid minor;
    legend('OAT used by JMOSS EKF','OAT derived from $T_{ic}$ and $K_t$','Location','NorthWest','Interpreter','latex');
    set(gcf,'Position',[0 0 1.6*windowSize windowSize]);
    title('\textbf{Temperature Estimation Performance, JMOSS V5 Demo}',...
    'Interpreter','latex','FontSize',fontSize+2)
    xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex','FontSize',fontSize); 
    ylabel('Estimated OAT, $T_{a}$, [K]','Interpreter','latex','FontSize',fontSize); 
    
    % Here let's take a look at the estimated variable recovery factor that
    % was produced by the JMOSS OAT estimator.
    subplot(3,1,2);
    plot(Mic,Kt,'b-','LineWidth',2);
    grid minor;
    xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex','FontSize',fontSize); 
    ylabel('Est. Recovery Factor, $K_{t}$','Interpreter','latex','FontSize',fontSize); 
    
    % Finally, let's look at the actual versus estimated Total Temprature
    % based on the JMOSS OAT estimator results.
    subplot(3,1,3);
    estimatedTotalTemp = estimatedOAT1.*(1 + 0.2*Kt.*Mic.^2);
    plot(Mic,Tic,'b-','LineWidth',2); hold on;
    plot(Mic,estimatedTotalTemp,'r-','LineWidth',2);
    legend('Actual $T_{ic}$','Estimated $T_{ic}$','Location','NorthWest','Interpreter','latex');
    grid minor;
    xlabel('Instrument Corrected Mach, $M_{ic}$','Interpreter','latex','FontSize',fontSize); 
    ylabel('Est. Total Temp., $T_{ic}$','Interpreter','latex','FontSize',fontSize);
    
end








