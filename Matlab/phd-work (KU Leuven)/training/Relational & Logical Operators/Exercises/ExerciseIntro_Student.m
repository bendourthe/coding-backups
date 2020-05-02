%% Basic operations using relational and logical operators
clear all
close all
clc
%% Part 1: Relational and logical operators
% Refresh your memory on relational operators (<, >, <=, >=, ==, ~=) using
% the following examples. Compare the results of the left and right
% syntax. Describe the syntax in words using commands %
a = 1; b = 2; c = 3; d = 1;
p1 = a == d;    p12 = eq(a,d);
p2 = c > b;     p22 = gt(c,b); 
p3 = c >= b;    p32 = ge(c,b); 
p4 = a <= c;    p42 = le(a,c); 
p5 = a ~= d;    p52 = ne(a,d);
p6 = a < d;     p62 = lt(a,d);

%% Part 2: Logical operators 
% More complex conditions can be represented by combining relational 
% operations using logical operators. The following logical operators are 
% available in MATLAB: ~, &, |, xor, &&, ||. Use the help function or check 
% your theory book. Compare the results of the left and right syntaxes. 
% Describe the syntax in words using commands %
A = [ 5, -3, 0, 0]; B = [2, 4, 0, 5];
z1 = A | B;                 z21 = or(A,B);
z2 = A & B;                 z22 = and(A,B);
z3 = ~A;                    z32 = not(A);
z4 = ~(A > 4); 			 
z5 = xor(A,B);
z6 = (c > b) && (a == d);
z7 = (c < b) && (a <= d);
z8 = (a < d) || (b > a);
z9 = any(A < 5 & B > 8);

% Note: The short circuit operators (&&, ||) will stop the evaluation of an 
% expression as soon as the results of the entire expression is known. Use 
% the short circuit operators  (&&, ||) when comparing single logical values 
% (scalars). Use the non-short circuit operators  (&, |) when comparing 
% arrays of logical values. Test the difference using the following
% commands:
Count = 0; Total = 1000;
(Count ~= 0) && (Total/Count > 80);             % Short-circuit
(Count ~= 0) & (Total/Count > 80);              % Non-Short-circuit
(rand(1,3) > rand(1,3)) & ([3 5 9] > [5 2 6]);

%% Part 3: Logical indexing
% This exercises show techniques of logical-indexing for matrix 
% manipulation.
x = [3 16 9 12 -1 0 -12 9 6 1];
y = [ 3 5 6 1 8 2 9 4 0 7];  
% Execute and interpret the results of the following commands:
a1 = x(x>6);
a2 = y(x<=4);
a3 = x(y<0);
a4 = (x>3)&(x<8);
a5 = x((x<2)|(x>=8));
a6 = y((x<2)|(x>=8));

% Create the following arrays using logical-indexing
% a) Set the positive values of x to zero in a new array x1;
x1 = x;
...
% b) Set values of x that are multiples of 3 to 3 (use rem) in a new array x2
x2 = x;
...
% c) Extract the values of x that are >10
x3 = x;
...

%% Part 4: combining relational and logical operators
% Create the following statements using relational and logical operators
C = randperm(10,10);% random permutations 
D = 1:10; 
% a) subtract from D taking 1 for C <= 7 and 0 otherwise 
...
% b) ones at positions where 2 <= C < 4 
...
% c) ones at positions where C == D or D == 3 
...
% d) 1 when ANY of the C elements are larger than 5 
...
% e) 1 when ALL D-elements are larger that 2 
...
% f) locate all nonzero elements of array C
...

%% Part 5: Leap year function
% A leap year is a year containing one additional day in order to keep the
% calendar year synchronized with the astronomical or seasonal year. To
% determine whether a year is a leap year or not you must check if the year
% is divisible by 400 or if the year is divisible by 4 and not by 100.
% Define a statement which checks if a year is leap year using logical
% operators and the rem function. Test the statement for different years.
year = 2013;
leapyear = [];...


