function [ASM, ekfResults] = JMOSSV31(TestPoints,varargin)
% Jurado-McGehee Online Self Survey Pitot-static Calibration Algorithm
% Authors: Juan Jurado and Clark McGehee
% Version 3.1: 21 March 2018
% Copyright 2018: Juan Jurado and Clark McGehee
%
% Usage:
% [ASM,ekfResults] = JMOSSV31(TestPoints)
%
% Required Inputs:
% TestPoints - k-dimensional structure containing source data with N
%              observations of required data.
%
%              TestPoints(k).Pt - Total pressure [psi]
%              TestPoints(k).Ps - Static presssure [psi]
%              TestPoints(k).Tic - Total temperature [K]
%              TestPoints(k).alphai - Indicated angle of attack [rad]
%              TestPoints(k).betai - Indicated angle of sideslip [rad]
%              TestPoints(k).Vn - North GPS speed [fps]
%              TestPoints(k).Ve - East GPS speed [fps]
%              TestPoints(k).Vd - Down GPS speed [fps]
%              TestPoints(k).hg - GPS altitude [ft MSL]
%              TestPoints(k).roll - Roll angle [rad]
%              TestPoints(k).pitch - Pitch angle [rad]
%              TestPoints(k).yaw - True heading angle [rad]
%              TestPoints(k).time - Absolute or relatve time vector [s]
%
% Optional Input Pairs:
% 'alpha' - Significance level for statistical inferences. Default 0.05,
%           which results in 95% inferences.
% 'maxKnots' - Maximum number of smoothing spline knots allowed. Default is
%              20. Supersonic data will have a minimum of 7.
% 'knots' - A P x 1 vector containing specific knot locations defined by
%           the user. This overrides the ASM automatic knot optimization 
%           algorithm and forces P knots at the specified Mic values.
% 'freezeStates' - true (default) or false boolean indicating whether or
%                  not to freeze the state estimates for wind, Kt, and dP0
%                  during the backwward pass. These states are usually
%                  stable on the backward pass without forcing them, and 
%                  are assumed to be constant. This setting forces them to 
%                  be constant in order to fully satisfy the assumption.
%
% Outputs:
% ASM - Akaike Spline Model object containing smooth Mach and dPp_Ps
%       results and model statistics.
%
%       ASM.mach = 1000 x 1 smooth vector spanning observed Mic domain
%       ASM.dPp_Ps = 1000 x 1 of corresponding dPp_Ps results 
%       ASM.PredictionBand = 1000 x 2 prediciton band for dPp_Ps
%       ASM.PredictionInterval = 1000 x 2 prediciton interval for dPp_Ps
%       ASM.ConfidenceBand = 1000 x 2 confidence band for dPp_Ps
%       ASM.ConfidenceInterval = 1000 x 2 confidence interval for dPp_Ps
%       ASM.maxPB = Maximum full width of Prediction Band
%       ASM.maxPI = Maximum full width of Prediction Interval
%       ASM.maxCB = Maximum full width of Confidence Band
%       ASM.maxCI = Maximum full width of Confidence Interval
%       ASM.splines = Mic values of any spline knots used in smoothing
%       ASM.statModel = MATLAB statistical model object for ASM
%
% ekfResults - Structure containing various EKF time histories for forward
%              and backward pass including estimates of 3D wind, Kt, and
%              non-standard pressure correction.
%
%              ekfResults(k).mach = N x 1 raw EKF observations of Mic for
%                                   k-th test point
%              ekfResults(k).dPp_Ps = N x 1 raw EKF estimates of dPp_Ps for
%                                     k-th test point
%              ekfResults(k).fwdPass = N x 6 array containing time
%                                      histories of the 6 EKF states 
%                                      (dPp_ps,Wn,We,Wd,Kt,dP0) for the
%                                      first pass on the k-th test point
%              ekfResults(k).bkdPass = N x 6 array containing time
%                                      histories of the 6 EKF states 
%                                      (dPp_ps,Wn,We,Wd,Kt,dP0) for the
%                                      second pass on the k-th test point
% 
% Requires:
% - MATLAB Version 2016B or later
% - Control System Toolbox
%
% Reference:
% Jurado, J.D., and McGehee, C.C., A Complete Online Algorithm for Air Data
%   System Calibration, AIAA Journal of Aircraft (Draft), March 2018.
%
% Copyright 2018: Juan Jurado and Clark McGehee

%% Feed inputs to driver function

    [ASM, ekfResults] = JMOSSDriver(TestPoints,varargin{:});    
end


