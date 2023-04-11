clear all; close all;
% Generate the initial floes, including the locations and radii
rng(100)
L = 36; % total number of the floes
r_min = 2*pi/100*4; % minimum radius
r_max = 2*pi/20*2; % maximum radius
radius = zeros(1,L);
l = 1;
% determining the floe radii
while l <= L
    radius_temp = randraw('pareto', [r_min,1],1); % power law random number
    if radius_temp <= r_max
        radius(l) = radius_temp;
        l = l + 1;
    end
end
% sort the array in a descend way, facilitating determining the locations
radius = sort(radius,'descend'); 
figure
plot(radius,'bo')
box on
set(gca,'fontsize',12)
title('Radius (non-dim)')
xlabel('Floe #')
% determining the floe locations
Location = zeros(2,L);
l = 1;
Location(:,l) = 2*pi*rand(2,1);
flag = 1;
for l = 2:L
    Diff_Threshold = radius(1:l-1) + radius(l);
    while flag > 0
        Location(:,l) = 2*pi*rand(2,1);
        Location_Diff_1a = abs(Location(1,l)*ones(1,l-1) - Location(1,1:l-1)); % distance in x
        Location_Diff_1b = 2*pi - Location_Diff_1a; % checking for periodic boundary
        Location_Diff_1 = min(Location_Diff_1a,Location_Diff_1b);
        
        Location_Diff_2a = abs(Location(2,l)*ones(1,l-1) - Location(2,1:l-1)); % distance in y
        Location_Diff_2b = 2*pi - Location_Diff_2a; % checking for periodic boundary
        Location_Diff_2 = min(Location_Diff_2a,Location_Diff_2b);
        

        Location_Diff = sqrt(Location_Diff_1.^2 + Location_Diff_2.^2);
        flag = sum(Location_Diff < Diff_Threshold-0.0);
    end
    flag = 1;
end
% Location(:,1) = [4,pi/2+0.5];
% Location(:,2) = [5,pi/2*3-0.5];

% floe thickness
thickness = 0.5*rand(1,L)*0 + 1; thickness = round(thickness * 100) / 100;

% plotting the floes
figure 
hold on
for l = 1:L
    th = 0:pi/50:2*pi;
    xunit = radius(l) * cos(th) + Location(1,l);
    yunit = radius(l) * sin(th) + Location(2,l);
    h = plot(xunit, yunit,'color',[0.2*radius(l),0.5,0.5]);
    text(Location(1,l),Location(2,l),num2str(thickness(l)));
end
xlim([0, 2*pi ])
ylim([0, 2*pi ])
box on
xlabel('x')
ylabel('y')
set(gca,'fontsize',12)
title('Floe locations (non-dim)')