function [x,fval] = optimN(a,b,c,d,e)
%     [X,FVAL] = optimS(A,B,C,D,E) solves an unconstrained
%   optimization problem in two dimension (old way of passing P1, P2).
%
%   Example
%  % 1) Smooth objective function
%     [x, fval ] = optimS(3,2,1,0,0)
%
%  % 2) Non-smooth objective function
%     [x, fval ] = optimS(3,2,1,5,0)
%
%  % 3) Stochastic objective function
%     [x, fval ] = optimS(3,2,1,0,2)
%
%  % 4) Non-smooth stochastic objective function
%     [x, fval ] = optimS(3,2,1,6,3)

% Copyright 2007 The MathWorks, Inc.

if nargin < 5
    a = 3; b = 2; c = 1; d = 6; e = 0;
end
range = [-10 10; -10 10];
x0 = [-5,9];     % Starting guess
options = optimset('LargeScale','off','OutputFcn',@myoutputfcn, ...
    'Display','off');
[x, fval] = fminunc(@objfun,x0,options);



    function f = objfun(x)
        % Nonlinear objective function
        f = a*x(:,1).^2 + b*x(:,1).*x(:,2) + ...
            c*x(:,2).^2 + d*abs(x(:,2) - x(:,1))+ e*randn;
    end

    function stop = myoutputfcn(x,optimvalues,state)
        stop = false;
        switch state
            case 'init'
                clf
                plotobjective(@(x) objfun(x),range);
            case 'iter'
                plot3(x(1),x(2),optimvalues.fval+2, ...
                    'ro','MarkerFaceColor','r','EraseMode','none', ...
                    'MarkerSize',6)
                pause(0.5);
            case 'done'
                plot3(x(1),x(2),optimvalues.fval+5, ...
                    'y*','MarkerFaceColor','y','EraseMode','none', ...
                    'MarkerSize',10)

        end
    end
end

