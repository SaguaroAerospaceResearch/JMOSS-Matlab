function [MODEL, EKF] = JMOSSV5(TestPoints,varargin)
% Jurado-McGehee Online Self Survey Pitot-static Calibration Algorithm
% Authors: Juan Jurado and Clark McGehee
% Version 5.0: 5 April 2019
% Copyright 2019: Juan Jurado and Clark McGehee
%
% Usage:
% [MODEL,EKF] = JMOSSV5(TestPoints)
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
%              TestPoints(k).hg - GPS altitude [ft HAE]
%              TestPoints(k).roll - Roll angle [rad]
%              TestPoints(k).pitch - Pitch angle [rad]
%              TestPoints(k).yaw - True heading angle [rad]
%              TestPoints(k).time - Absolute or relatve time vector [s]
%
% Note: 
% You may provide OAT measurements (instead of Total Temperature) using:
%
%              TestPoints(k).Ta - Ambient temperature [K]
%
% but be sure to remove the Tic field, otherwise Ta will be ignored.
%
% Optional Input Pairs:
% 'alpha' - Significance level for statistical inferences. Default 0.05,
%           which results in 95% inferences.
% 'rollLimit' - The Angle of Bank (AOB) limit for excluding turn data from
%               dPp_Ps modeling/smoothing. Default is 10 [deg] AOB.
%
%
% Outputs:
% MODEL - Gaussian Regression Process (GRP) object containing smoothed
%         instrument corrected Mach number (Mic) and corresponding static 
%         position error ratio (deltaPp_Ps).
%
%       MODEL.mach = 1000 x 1 smooth vector spanning observed Mic domain
%
%       MODEL.dPp_Ps = 1000 x 1 of corresponding dPp_Ps results 
%
%       MODEL.PredictionBand = 1000 x 2 prediction band for dPp_Ps
%
%       MODEL.maxPB = Maximum full width of Prediction Band
%
%       MODEL.fullModel = Full MATLAB GRP object
%
%       MODEL.predictionFun = The final dPp_Ps function as a function of
%       input machIC. This function outputs the best estimate for dPp_Ps
%       for any given input machIC number(s).
%
%       MODEL.Curves = Structure containing 1000 x 1 smooth vectors for 
%       deltaH, deltaV, deltaM, and their corresponding 1000 x 2 prediction
%       bands. These are correction curves and assume sea level.
% 
% EKF - Structure containing various EKF time histories for forward
%              and backward pass including estimates of 3D wind, r, and
%              pressure linearization point.
%
%              ekfResults(k).mach = N x 1 raw EKF observations of Mic for
%              k-th test point
%
%              ekfResults(k).machPC = N x 1 raw EKF observations of pressure
%              corrected Mach number for k-th test point
%
%              ekfResults(k).dPp_Ps = N x 1 raw EKF estimates of dPp_Ps for
%              k-th test point
%
%              ekfResults(k).fwdPass = N x 5 array containing time
%              histories of the 5 EKF states (dPp_ps,Wn,We,Wd,dP0) for the
%              first pass on the k-th test point
%
%              ekfResults(k).bkdPass = N x 5 array containing time
%              histories of the 5 EKF states (dPp_ps,Wn,We,Wd,dP0) for the
%              econd pass on the k-th test point
%
%              ekfResults(k).deweightedPts = N x 1 boolean indexing vector
%              indicating which samples were included (true) or excluded
%              (false) from the GRP smoothing process based on rollLimit.
%
%              ekfResults(k).subSamp = N x 1 vector containing the indices
%              of the samples that were used in the EKF process after
%              decimating the input data to remove repeated measurements.
%
%              ekfResults(k).Ta = N x 1 vector containing the JMOSS
%              estimate of OAT, for the k-th test point, based on Tic. 
%              This field is only created if the user provided a Tic 
%              (Total Temperature) input. 
%
%              ekfResults(k).Kt = N x 1 vector containing the JMOSS
%              estimate of the temperature recovery factor, Kt, for the
%              k-th test point. This field is only created if the user
%              provided a Tic (Total Temperature) input.
%
% Requires:
% - MATLAB Version 2015B or later
% - Optimization Toolbox 
% - Statistics and Machine Learning Toolbox (for GRP smoothing; can still
%   run EKF process without it).
%
% Reference:
% Jurado, J.D., and McGehee, C.C., A Complete Online Algorithm for Air Data
%   System Calibration, AIAA Journal of Aircraft, 29 October 2018.
%
% Copyright 2019: Juan Jurado and Clark McGehee

%% Feed inputs to driver function

    [MODEL, EKF] = JMOSSDriver(TestPoints,varargin{:});    
end


