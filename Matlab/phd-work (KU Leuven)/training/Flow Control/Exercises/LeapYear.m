%% Leap year check
% A leap year is a year containing one additional day in order to keep the
% calendar year synchronized with the astronomical or seasonal year. To
% determine whether a year is a leap year or not you must check if the year
% is divisible by 400 or if the year is divisible by 4 and not by 100.
% Define a statement which checks if a year is leap year using logical
% operators and the rem function. Test the statement for different years.
year = 2000;
if (rem(year,4) == 0)
    if (rem(year,100) == 0)
        if (rem(year,400) == 0)
            LeapYear = true;
        else
            LeapYear = false;
        end
    else
        LeapYear = true;
    end
else
    LeapYear = false;
end
% Print leap year output in command window
[year LeapYear]

% Another option
LeapYear2 = (rem(jaar,4) == 0) && (rem(year,100) ~= 0 || rem(year,400) == 0);
% Print leap year output in command window
[year Leapyear2]

%% Days in a month
% Write a switch block to check how many days are in a month. In a leap 
% year Februari contains 29 days
year = 2000; % Year
month = 4;   % Month
switch month
    case {4,6,9,11}
        nDays = 30;
    case 2
        if (rem(year,4) == 0) && (rem(year,100) ~= 0 || rem(year,400) == 0)
            % It is a leap year
            nDays = 29;
        else
            nDays = 28;
        end
    otherwise
        nDays = 31;
end
[year month nDays]  % Display result in command window