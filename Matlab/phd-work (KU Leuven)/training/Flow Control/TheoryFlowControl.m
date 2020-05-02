%% Theory examples
% The purpose of this part of the course is to make you familiar with
% - IF statements
% - SWITCH statements
% - FOR loops
% - WHILE loops

%% If-elseif-else statements
% An if statement evaluates an expression, and executes a group of
% statements when the expression is true. An evaluated expression is true
% when the result is nonempty and contains all nonzero elements (logical or
% real numeric). Otherwise the expression is false. The expressions
% containing often relational and logical operators.
z = [3,5];
if z(1) > 2
    z1_gt_2 = 1;
end

% If the statement is false, nothing will happen. If you do want something
% to happen, use else. The else statements are executed only when previous
% expressions in the IF block are false.
if z(1) > 4;
    z1_gt_4 = 1;
else
    z1_gt_4 = 0;
end

% To check for more expressions, use elseif. It is possible to include
% multiple elseif statements in an if block. When one expression is
% verified, the statements associated with it are executed and the if-block
% does not continue verifying the rest of the expressions.
if z(1) > 3;
    z1_gt_3 = 1;
elseif z(2) > 4;
    z2_gt_4 = 1;
else
    z1_gt_3 = 0;
    z2_gt_4 = 0;
end

% More examples using the if-construct.
B = -10;
if B < 0
    disp('B is less than 0');
end

% Here is a code that checks the age of a person:
age = 17;
if age <= 18;
    disp('Person is 18 or younger');
else
    disp('Person is older than 18');
end

% If you want to test for more age ranges; e.g. 18 or younger,
% 19 to 21, and 22 and over you can include an elseif construcion.
age = 17;
if age <= 18;
    disp('Person is 18 or younger');
elseif age < 22;
    disp('Person is between 19 and 21');
else
    disp('Person is older than 21');
end

% Find the lowest number using the if structure
A = 237; B = 53; C = 94;
if (A < B) && (A < C)
    Result = A;
elseif (B < A) && (B < C)
    Result = B;
else
    Result = C;
end

%% Switch-case-otherwise statement
% A switch block conditionally executes one set of statements from several
% choices. Each choice is a case.

% Numeric expressions
time = 0:0.1:2*pi;
signal = [sin(time);cos(time);2*sin(time)];
n_signal = 1;   % signal number
switch n_signal
    case 1
        out_sig = signal(1,:);
    case 2
        out_sig = signal(2,:);
    case 3
        out_sig = signal(3,:);
    otherwise
        out_sig = [];
end

% String expressions
X = rand(30,4);
plottype = 'bar';
switch plottype
    case '2Dline'
        plot(X)
        title('2D line graph')
    case 'box'
        boxplot(X)
        title('Box plot')
    case 'bar'
        bar(X)
        title('Bar graph')
    otherwise
        warning('Unexpected plot type');
end

% Cell array case expression
month = 10;
switch month
    case {12, 1, 2}
        disp([num2str(month) ' is a winter month'])
    case {3, 4, 5}
        disp([num2str(month) ' is a spring month'])
    case {6, 7 , 8}
        disp([num2str(month) ' is a summer month'])
    case {9, 10, 11}
        disp([num2str(month) ' is a fall month'])
    otherwise
        disp([num2str(month) ' is not a valid month'])
end

Month = 'Jan';
switch Month
    case {'Dec','Jan','Feb'}
        disp([num2str(Month) ' is a winter month'])
    case {'Mar','Apr','May'}
        disp([num2str(Month) ' is a spring month'])
    case {'Jun','Jul','Aug'}
        disp([num2str(Month) ' is a summer month'])
    case {'Sep','Oct','Nov'}
        disp([num2str(Month) ' is a autumn month'])
    otherwise
        disp([num2str(Month) ' is not a valid month'])
end

%% For loop
% In this for loop i is increased by 1 every step until 10
for i = 1:10
    step = i;
end

% In this for-loop the index is increased by 1 for every step so
% index = 3, 4, 5, 6
for index = 3:6
    H = index*2;
    disp(H)
end

% For statements can also be very useful if you want to do the same
% calculations over different files.
for k = 1:5
    file = ['subject' num2str(k) '.mat'];
end

% Fibonacci
% Initialize the first two values
f(1) = 0; f(2) = 1;
% Create the first 30 Fibonacci numbers
for i = 3:30
    % Perform the sum of terms accordingly
    f(i) = f(i-1) + f(i-2);
end

% Row: Tn = sum_j=1 till n (-1)^j-1*1/j
n = 10;
sign = -1;
a = 0;
for j = 1:n
    sign = -sign;
    a = a + sign.*1./j;
end
Tn = a;

% Instead of increasing the increment every step by 1 you can define the
% stepsize. This example increases the index value by 3 for every step.
for index = 3:3:12
    H = index*2;
    disp(H)
end


y = [1,5,4,5,6; 7,6,2,1,4; 5,8,1,3,9; 8,5,4,2,7];% Example matrix
av_coly = mean(y,1);    % average column value of matrix y
% If you want compute the average value of a column for every second
% column. Note that the index j is used to define the subscript indices of
% the variable av_coly2
for j = 1:2:4
    av_coly2(j) = mean(y(:,j));
end
% For large samples, you first want to know the matrix size
nrow = size(y,1);
for i = 1:2:nrow
    av_rowy2(i) = mean(y(i,:));
end

% Figure with different sine functions with an increment size of 3
x = 2*pi*(0:0.01:1);
figure, hold on
colors = {'r','g','b','y','m','c','k'}; % line color
i = 1;                  % start index color
for a = 1:3:9
    plot(x, sin(a*x),'color',colors{i})
    i = i + 1;          % updat index color
end

% Figure with different line colors using a cell construct to define the
% steps of the for-loop
figure, hold on
for colors = {'r','g','b','y','m','c','k'}; % line color
    b = rand([100,1]);
    plot(b, 'Color',colors{1})
end

% The steps of this for-loop are defined by the matrix x
for x = [4 7 2 -4 10]
    y = x^2 + x;
    disp('Current value of y:')
    disp(y)
end

% Nested for loop
% Define the dimensions of a matrix
nr = 4; % number of rows
nc = 3; % number of columns
p = 1;  % start value
A=zeros(nr,nc);
for r=1:nr
    for c=1:nc
        A(r,c) = p;
        p = p+1;
    end
end

%% Avoid for loops
y=randn(500,500);

% Solution with for loop
tic
[row,col,v] = find(y>0);
y(row,col)=0;
for i = 1:size(row);
    y(row(i),col(i))=0;
end
toc % Elapsed time is 13.703791 seconds.

% Solution with logical indexing
tic
y(y>0)=0;
toc % Elapsed time is 0.005004 seconds.

%% Combining for loops with if and/or switch statements
C = [7,12,9,18,14,6];
for i = 1:length(C)
    if C(i) <= 10
        D(i) = 1;
    elseif C(i) <= 15;
        D(i) = 2;
    else
        D(i) = 3;
    end % if statement
end % for loop

%% While loops
% A simple example of a while loop is to find the first integer n for which
% its factorial is a three-digit number.
n = 1;          % initial n value
nFactorial = 1; % initial factorial number
while nFactorial < 1e3
    nFactorial = nFactorial * n;
    disp('Current value of nFactorial:')
    disp(nFactorial)
    n = n + 1; 	% update n value
end
n_Factorial = n - 1; % Factorial of n

% Make a m-script which divide an arbitrary number x by 2 till the new
% obtained number is smaller than 2. The output of the script is a number p
% that represents how often the number x is divided by 2 and the value q
% that is left when the number x is divided by 2, p times.
x = 49;
p = 0;
while x > 2
    p = p + 1;  % update p
    q(p) = x/2;
    x = q(p);
end

% Recurrent formula
x = 0:0.001:2;
y = zeros(size(x));
y(1) = 1;
i = 1;
while(i < length(x))
    y(i+1) = y(i) + h*(x(i)-abs(y(i)));
    i = i + 1;
end
plot(x,y)

% Nested while-loop
nr = 4; % number of rows
nc = 3; % number of columns
i = 1;  % set start integer i
j = 1;  % set start integer j
while i<=nr
    while j<=nc
        H(i,j) = i+j;
        j = j+1;  % update j value
    end    % end inner loop
    j = 1;         % set j back to 1
    i = i + 1;     % update i value
end         % end outer loop

%% Continue
% The continue statement passes control to the next iteration of the for or
% while loop in which it appears, skipping any remaining statements in the
% body of the loop.
B = [3, 5, 7, 0, 20];
ProductB = 1;
for i = 1:length(B)
    if B(i) == 0 % When B(i) equals to zero, it skips the statement ProductB = ProductB * B(i)
        continue;
    end
    ProductB = ProductB * B(i);
end

% However, the continue statement is rarely needed and its use should be
% avoided. For example the previous code example could be written as
B = [3,5, 7, 0, 20];
for i = 1:length(B)
    if B(i) ~= 0
        ProductB = ProductB * B(i);
    end
end

%% Break
% The break statement terminates the execution of a for or while loop
% before all predefined conditions are met. The remaining statements in the
% loop that appear after the break statement are not executed.
b = 0;
while 1
    b = b + 1;
    if b > 5
        break;
    end
end
