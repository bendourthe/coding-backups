classdef blip < handle
    %BLIP Class for representing a radar blip
    %   A blip knows its angle, signal and range

    % Copyright The MathWorks, Inc. 2008, 2010

    properties
        AoA
        range
        signal
    end

    methods
        function obj = blip(AoA, range, signal)
            if nargin == 3
                obj.AoA    = AoA;
                obj.range  = range;
                obj.signal = signal;
            end
        end
        
        function identify(obj)
            disp(['I am a friend, please don''t shoot!' ...
                'I am arriving at ' num2str(obj.AoA)]);
        end
    end
end
