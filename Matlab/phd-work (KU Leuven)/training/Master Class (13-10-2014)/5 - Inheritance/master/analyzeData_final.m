%% Analyze Sensor Design
% We explore a sensor design for a far-field staring array.
% The array is designed to detect the angle of arrival (AoA)
% for a distant object emitting a signal of some frequency.
%
% The system is tested with multiple targets.  Different
% sensor configurations can be explored to investigate the
% accuracy of AoA detection, as well as the ability of the
% sensor to resolve multiple targets.

% Copyright The MathWorks, Inc. 2008, 2010

%% Clear Workspace
clear all
close all

%% Define Radio Beacon
signal.spdOfLight  =    3e8 ;
signal.freq        =   12e7 ;
signal.wavelength  = signal.spdOfLight/signal.freq  ;
signal.amp         =      1 ;

%% Define Aircraft
balloon{1}          = blip(-10, 1e8, signal) ;
balloon{2}          = blip( 20, 1e8, signal) ;
balloon{3}          = movingBlip( 5, 25, 1e8, signal) ;

%% Define Sensor Array
staringArray      = sensor ;

%% Compute Angles of Arrival
arrivalAngles      = staringArray.computeAoA(balloon) ;
disp(arrivalAngles)

%% Move the target
for count = 1:15
    balloon{3}.move() ;
    arrivalAngles  = staringArray.computeAoA(balloon) ;
    pause(0.5) ;
end