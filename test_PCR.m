% Createdby Ying Du
% All rights reserved

close all
clear all
clc


% this code simulates the temperature distribution along the tube

v = 40; % unit: uL/sec
a = 0.5; % width of the channel, unit: mm
b = 0.5; % height of the channel, unit: mm
L = 16; % length of the channel, unit: mm

T_liquid = 95; % initial temperature of the liquid, unit: °C
T_wall = 25; % temperature of the surrounding wall, unit: °C

alpha = 0.143; % thermal diffusivity of liquid, unit: mm/s

%% user should NOT change from here

T_delta = T_liquid - T_wall; % initial temeprature difference bettwen the liquid and the wall
velocity = v/a/b; % velocity of the liquid in the channel, mm/s

c = sqrt(alpha);

% numerical simulation parameters
M_limit = 100;
N_limit = 100;


[M, N] = meshgrid(1:M_limit, 1:N_limit);

A_mn = 4 * T_delta / pi^2 ./ (M .* N) .* (1 - (-1).^M) .* (1 - (-1).^N);
lambda_mn = c * pi *sqrt(M.^2/a^2 + N.^2/b^2);

x = linspace(0, a, 100);
y = linspace(0, b, 100);
l = linspace(0, L, 256);
t = l/velocity;


% claculate cross-section on x-z plane
[L_xz, X_xz] = meshgrid(l, x);

T_xz = zeros(size(L_xz));

temp_y = b/2; % we choose the cross section at y = b/2;

T_xz = T_xz + T_wall;
for m = 1:M_limit
    for n = 1:N_limit
        T_xz = T_xz + A_mn(n, m) .* sin(m*pi/a*X_xz) .* sin(n*pi/b*temp_y) .* exp(-lambda_mn(n,m)^2*L_xz/velocity);
    end
end
 
% claculate cross-section on y-z plane
[L_yz, Y_yz] = meshgrid(l, y);

T_yz = zeros(size(L_yz));

temp_x = a/2; % we choose the cross section at x = a/2;

T_yz = T_yz + T_wall;
for m = 1:M_limit
    for n = 1:N_limit
        T_yz = T_yz + A_mn(n, m) .* sin(m*pi/a*temp_x) .* sin(n*pi/b*Y_yz) .* exp(-lambda_mn(n,m)^2*L_yz/velocity);
    end
end
 
figure
subplot(2,1,1)
imagesc(l, x, T_xz)
set(gca, 'fontsize', 16)
xlabel('Length (mm)', 'fontsize', 16)
ylabel('Width (mm)', 'fontsize', 16)
title(['Width-length plane, ' num2str(v) ' \muL/sec'], 'fontsize', 16)
caxis([min([T_liquid T_wall]) max([T_liquid T_wall])]);
colormap jet
colorbar
grid on

subplot(2,1,2)
imagesc(l, y, T_yz)
set(gca, 'fontsize', 16)
xlabel('Length (mm)', 'fontsize', 16)
ylabel('Heigt (mm)', 'fontsize', 16)
title(['Height-length plane, ' num2str(v) ' \muL/sec'], 'fontsize', 16)
caxis([min([T_liquid T_wall]) max([T_liquid T_wall])]);
colorbar
grid on


% plot
T_center = T_xz(round(length(x)/2), :);
figure
plot(l, T_center, 'linewidth', 3);
set(gca, 'fontsize', 16)
xlabel('Length (mm)', 'fontsize', 16)
ylabel('Temperature (°C)', 'fontsize', 16)
title(['Temperature vs. Channel length, ' num2str(v) ' \muL/sec'], 'fontsize', 16)
grid on

% plot
figure
plot(t, T_center, 'linewidth', 3);
set(gca, 'fontsize', 16)
xlabel('Time (sec)', 'fontsize', 16)
ylabel('Temperature (°C)', 'fontsize', 16)
title(['Temperature vs. time, ' num2str(v) ' \muL/sec'], 'fontsize', 16)
grid on