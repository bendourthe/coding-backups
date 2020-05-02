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
signal.wavelength  = signal.spdOfLight/signal.freq ;
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

%% Create a listener to an Event
% listen to an event issued with the notify command
listen             = addlistener(balloon{1},'blipMoved', ...
    @(varargin)disp('Notified of an event'));

%% Trigger notify for event
% moveBlip method calls notify
balloon{1}.moveBlip(-60);

%% Create a listener to an Observable Property
% The blip class has SetObservable property on AoA
listen2            = addlistener(balloon{2},'AoA','PostSet', ...
    @(src,evnt)disp('Observable Property was Changed'));

%% Observable Property automatically triggers event
balloon{2}.AoA = -45;

%% Register Blips with Sensor Array
staringArray       = staringArray.registerTarget(balloon);

%% Start the movingBlip moving
balloon{3}.start() ;

%% Stop the blip
balloon{3}.stop() ;

%% Clean up
close all;
clear all;
delete(timerfindall);
clc;