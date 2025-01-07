clc
clear all
close all
%Parameters
b = 1 %width
h = 2% height
E = 2e11 % Young's Modulus (Pa) for steel
I = (b*h.^3)/12; %Moment of Inertia (m^4)
w = 200000000; %kN
l = 10; %meters
x = linspace(0,l,1000); %1000 points along the beam
V = (-(w*x)+(w*l)/(2));
M = (-(w*x.^2)/2 + (w*l*x)/2);

% Calculate second derivative of displacement (u'')
d2u_dx2 = M / (E * I);

%Integrate to find slope (u')
slope = cumtrapz(x, d2u_dx2);

%Integrate again to find displacement (u)
u = cumtrapz(x,slope);

% Apply boundary conditions (adjust displacement to set u(0) = 0 and u(L) = 0)
u = u - u(1); %Ensure u(0)=0
u = u - u(end)*x/l %ensure u(L)=0

%plot u displacement
figure;
plot(x,u,'k-');
grid on;
xlabel('Position along the beam (m)');
ylabel('Displacement (m)');
title('Displacement of Simply Supported Beam');
[minValue, minIndex] = min(u)
xMin = x(minIndex);
hold on;
plot(xMin, minValue, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red circle
text(xMin, minValue, sprintf('Min: (%.2f, %.2f)', xMin, minValue), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'Color', 'red'); % Add text annotation
hold off;

%plot shear force
figure
plot (x,V)
grid on;
xlabel ('x(m)')
ylabel ('Shear Force (kN)')
title ('Shear force diagram')
hold on;
yline (0, 'k-', 'LineWidth', 0.5)
[maxValue, maxIndex] = max(V)
xMax = x(maxIndex);
plot(xMax, maxValue, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red circle
text(xMax, maxValue, sprintf('Max: (%.2f, %.2f)', xMax, maxValue), ...
    'VerticalAlignment', 'cap', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'Color', 'red'); % Add text annotation
[minValue, minIndex] = min(V)
xMin = x(minIndex);
plot(xMin, minValue, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red circle
text(xMin, minValue, sprintf('Min: (%.2f, %.2f)', xMin, minValue), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'Color', 'red'); % Add text annotation
hold off;

%plot bending moment
figure
plot(x,M)
grid on;
xlabel ('x(m)')
ylabel ('Bending Moment (kN-m)')
title ('Bending Moment diagram')
hold on;
[maxValue, maxIndex] = max(M)
xMax = x(maxIndex);
plot(xMax, maxValue, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Red circle
text(xMax, maxValue, sprintf('Max: (%.2f, %.2f)', xMax, maxValue), ...
    'VerticalAlignment', 'cap', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'Color', 'red'); % Add text annotation
hold off;
