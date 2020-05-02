classdef blip
    %BLIP Class for representing a target
    %   A blip class knows its angle, signal and range

    % Copyright The MathWorks, Inc. 2008, 2010

    properties
        AoA
        range
        signal
    end

    methods
        function obj = blip(AoA, range, signal)
            if nargin == 3
                obj.AoA    = AoA   ;
                obj.range  = range ;
                obj.signal = signal;
            end
        end
        function identify(obj)
            disp(['I am a friend, don''t shoot. ' ...
                'I am arriving at ' num2str(obj.AoA)]);
        end
    end
end
