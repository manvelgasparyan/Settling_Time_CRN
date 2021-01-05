%--------------------------------------------------------------------------
%The stoichiometric matrix
S = [-2  2  0  0;
      1 -1  0  0;
      7 -5 -1  1;
      0  0  5 -5];
%--------------------------------------------------------------------------
%The lengh of the time interval
T = 1.55*10^(-3);
%--------------------------------------------------------------------------
%Number of time points
mu = 1000;
%--------------------------------------------------------------------------
%Initial concentrations
x_0 = [5; 0; 4; 0];
%--------------------------------------------------------------------------
%Rate constants
k = [127; 40.5; 26; 10.6];
%--------------------------------------------------------------------------
%Michaelis constants
K = [100; 50; 20; 45; 50];
%--------------------------------------------------------------------------
%Species' concentrations
[t,x] = Concentrations (k, K, S, mu, T, x_0);
%--------------------------------------------------------------------------
%Number of species
[s, ~] = size(x);
%--------------------------------------------------------------------------
figure(1)
%delta-bound function
p = patch([0 0 1.6*10^(-3) 1.6*10^(-3)],[4.997 5.001 5.001 4.997],'g');
%---------------------------------
hold on
%Concentration
plot(t, x(1,:), 'color', 'b', 'linewidth', 1.5)
%----------------------------------
hold on
%Steady State
x_star = [0 1.6*10^(-3)];
y_star = [4.999 4.999];
line(x_star, y_star, 'Color' , 'r' , 'LineStyle' , '--' , 'LineWidth', 1.5)
%----------------------------------
hold on
%Settling time
x_set = [1.054*10^(-3) 1.054*10^(-3)];
y_set = [0 5.01];
line(x_set, y_set, 'Color' , 'black' , 'LineStyle' , '--' , 'LineWidth', 1.5)
%----------------------------------
hold on
plot(1.054*10^(-3), 4.997,'r*')
%----------------------------------
hold off
grid on
grid minor     
caption3 = sprintf('Concentration');
set(gca,'FontSize',12)
xlabel({'Time','(in seconds)'})
ylabel('Concentration')
title(caption3)
ylim([4.9 5.05])
%--------------------------------------------------------------------------