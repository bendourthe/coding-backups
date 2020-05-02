classdef blip
    %BLIP Definition of all radar blips
    %   This is a radar blip class for sensor modelling

    % Copyright The MathWorks, Inc. 2008, 2010

    properties
        AoA
        range
        signal
    end

    methods
        function obj = target(AoA, range, signal)
            if nargin == 3
                obj.AoA    = AoA    ;
                obj.range  = range  ;
                obj.signal = signal ;
            end
        end
        
        function identify(obj)
            disp(['I am a friend, please don''t shoot. ', ...
                'I am arriving at ' num2str(obj.range)]) ;
        end
    end
end
