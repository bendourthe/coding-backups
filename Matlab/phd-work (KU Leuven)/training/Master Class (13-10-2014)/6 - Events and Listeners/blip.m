classdef blip < handle
    %AIRCRAFTCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetObservable = true)
        AoA
    end
    
    properties
        range
        signal
    end
    
    events
        blipMoved
    end
    
    methods
        function obj = blip(AoA, range, signal)
            if nargin == 3
                obj.AoA    = AoA    ;
                obj.range  = range  ;
                obj.signal = signal ;
            end
        end
        
        function obj = moveBlip(obj,AoA)
            obj.AoA    = AoA    ;
            notify(obj, 'blipMoved');
        end
    end
    
end
