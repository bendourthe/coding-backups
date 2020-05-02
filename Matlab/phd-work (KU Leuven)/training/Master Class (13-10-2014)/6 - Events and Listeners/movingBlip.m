classdef movingBlip < blip
    %MOVINGBLIP Summary of this class goes here
    %   Detailed explanation goes here

    % Copyright The MathWorks, Inc. 2008, 2010

    properties
      timerObj
      deltaAoA  = 5;
    end

    methods
        %% Constructor
        function obj = movingBlip(deltaAoA, varargin)
            % assign the superclass portion
            obj = obj@blip(varargin{:}) ;
            
            if nargin >= 1
                % assign the movingBlip's unique property
                obj.deltaAoA = deltaAoA ;
            end
            
            % Configure a timer object to move blip
            obj.timerObj = timer('TimerFcn'     , @obj.timerMove , ...
                                 'ExecutionMode', 'fixedrate'    , ...
                                 'Period'       ,  0.1           );

        end

        %% Move function to explicitly move blip
        function move(obj)
            % add change in AoA
            obj.AoA =  obj.AoA + obj.deltaAoA ;

            % Keep blip in in bounds by flipping direction after having
            % moved too far
            if abs(obj.AoA) > 65
                obj.deltaAoA = -obj.deltaAoA ;
            end

        end
        
        %% Timer Start and Stop Functions
        function start(obj)
            start(obj.timerObj);
        end
        
        function stop(obj)
            stop(obj.timerObj);
        end
        
        %% Destructor to insure timer is deleted with object
        function delete(obj)
            stop(  obj.timerObj);
            delete(obj.timerObj);
        end
        
        function timerMove(obj, varargin)
            obj.move();
        end
    end
end
