%% Dice
% Dice are used for generating random numbers and are often used in 
% gambling games. A traditional dice is a cube containing numbers 1-6.
% The following statement allows you to create random numbers (integers)
% from 1 - 6.
Dice = randi(6,1);

% Create a program that count the number of rolls before you roll 6. 
% (Tip: while statement)
nRoll = 0;
Dice = 0;
while Dice ~= 6
    nRoll = nRoll + 1;
    Dice = randi(6,1);
end

% Adapt the program such you can count the number of rolls needed to roll 5 
% for 10 subjects.
nRoll2 = zeros(1,10);
Dice2 = zeros(1,10);
for i = 1:lenght(nRoll2)
    while Dice2(i) ~= 5
        nRoll2(i) = nRoll(i) + 1;
        Dice2(i) = randi(6,1);
    end
end