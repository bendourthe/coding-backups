classdef movingBlip < blip
    %MOVINGBLIP Summary of this class goes here
    %   Detailed explanation goes here

    % Copyright The MathWorks, Inc. 2008, 2010

    properties
        deltaAoA
    end

    methods
        function obj = movingBlip(deltaAoA, varargin)
            % assign the superclass portion
            obj = obj@blip(varargin{:}) ;
            
            if nargin >= 1
                % assign the movingBlip's unique property
                obj.deltaAoA = deltaAoA ;
            end
        end

        function move(obj)
            % add change in AoA
            obj.AoA =  obj.AoA + obj.deltaAoA ;

            % Keep blip in in bounds by flipping direction after having
            % moved too far
            if abs(obj.AoA) > 65
                obj.deltaAoA = -obj.deltaAoA ;
            end

        end

    end
end
