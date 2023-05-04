% 2D shallow water equation
K_max = 2; % the range of Fourier modes is [-K_max, K_max]^2
k = zeros(2, (2 * K_max + 1) * (2 * K_max + 1)); % Total number of Fourier wavenumbers

% arranging Fourier wavenumbers
% arranging in such a way that the complex conjugates modes are next to
% each other, namely (-k1,-k2) will be next to (k1,k2). This will
% facilitate data assimilation and allows it to be a nearly block diagonal
% matrix where each block is 2 by 2.

m = 1;
for i = - K_max : K_max
    if i < 0
        for j = - K_max: i
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    else
        for j = - K_max : i - 1
            k(1, m) = i;
            k(2, m) = j;
            m = m + 2;
        end
    end
end
k(:, 2: 2: end - 1) = - k(:, 1 : 2 : end - 2);

% showing the grids of Fourier wavenumbers
figure
hold on
plot(k(1,1:2:end-2), k(2,1:2:end-2), 'ro', 'linewidth',4);
plot(k(1,2:2:end-1), k(2,2:2:end-1), 'go', 'linewidth',4);
plot(k(1,end), k(2,end), 'ko', 'linewidth',4);
box on
grid on
set(gca, 'fontsize', 12)
title('Fourier wavenumbers', 'fontsize', 14)
xlabel('k_1')
ylabel('k_2')

% Total number of Fourier modes is 3 times as many as that of the Fourier
% wavenumbers because each Fourier wavenumber has three modes: one GB mode
% that is incompressible and two gravity modes that are compressible
kk_end = [0;0];
kk = [k,k,kk_end,k,kk_end]; % total number of the waves, 2k gravity, k-1 GB and two background modes
kk(:,length(k(1,:))) = []; % putting the two +/-(0,0) modes of the gravity waves together
epsilon = 0.2; % Rossby number
% dispersion relationship for the two gravity modes, p and m stand for plus and minus
% the dispersion relationship for the GB mode is omega = 0
omegak_p = 1/epsilon * sqrt(k(1,:).^2 + k(2,:).^2 + 1); 
omegak_m = - 1/epsilon * sqrt(k(1,:).^2 + k(2,:).^2 + 1);
omegak = reshape([omegak_p; omegak_m],1,[]);
% eigenvectors of GB and gravity modes
% the last column of rk2 and rk3 are the (0,0) gravity modes, which need to
% deal with in a different way
% rk1: GB; rk2: gravity +; rk3: gravity -.
rk1 = [1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (-1i * k(2,:));
    1./sqrt(k(1,:).^2 + k(2,:).^2+1) .* (1i * k(1,:))];
rk1(:,end) = [1;0]; rk1_end = [0;1]; rk1 = [rk1,rk1_end];
rk2 = [1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)./sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (1i * k(2,:) + k(1,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1));
    1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)./sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (-1i * k(1,:) + k(2,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1))];
rk2(:,end) = [1i;1]/sqrt(2);
rk3 = -[1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)./sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (1i * k(2,:) - k(1,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1));
    1./sqrt(k(1,:).^2 + k(2,:).^2)/sqrt(2)./sqrt(k(1,:).^2 + k(2,:).^2 + 1) .* (-1i * k(1,:) - k(2,:) .* sqrt(k(1,:).^2 + k(2,:).^2 + 1))];
rk3(:,end) = [-1i;1]/sqrt(2);

rk = zeros(size(kk));
rk(:, 1:2:length(k(1,:))*2-1) = rk2;
rk(:, 2:2:length(k(1,:))*2) = rk3;
rk(:, length(k(1,:))*2+1:end) = rk1;


% stochastic systems for each Fourier mode
% written in the vector form
% dU = (a1 U + a0 ) dt + b1 dW

rng(500)
T = 10; % total time
dt = 0.001; % time step
N = round(T/dt);
Dim_U = length(kk(1,:)); % dimension of the system
Dim_Ug = length(k(1,:)); Dim_UB = Dim_Ug - 1;
u_hat = zeros(Dim_U,N); % define all the Fourier modes
d_B = 0.5; % damping of the GB modes
d_g = 0.5; % damping of the gravity modes
sigma_B = 0.4; % noise of the GB modes
sigma_g = 0.4; % noise of the gravity modes
f_amp = 1.6; % large-scale forcing amplitude
f_phase = 1.5; % large-scale forcing period
f_x_b = 0.5; % large-scale forcing background in x direction
f_y_b = 0;% large-scale forcing background in y direction
% b1: noise coefficient; a1: damping and phase 
b1 = zeros(Dim_U, Dim_U);
a1 = - diag([d_g * ones(1,Dim_Ug * 2), d_g * ones(1,Dim_UB), d_B, d_B]) + 1i * diag([omegak, zeros(1,Dim_UB+2)]);
L_u = a1;
for i = 1: 2: 2 * Dim_Ug + Dim_UB
    if i <= 2 * Dim_Ug
        b1(i,i) = 1 / sqrt(2) * sigma_g;
        b1(i+1, i+1) = -1i / sqrt(2) * sigma_g;
        b1(i, i+1) = 1i / sqrt(2) * sigma_g;
        b1(i+1, i) = 1 / sqrt(2) * sigma_g;
    else
        b1(i,i) = 1 / sqrt(2) * sigma_B;
        b1(i+1, i+1) = -1i / sqrt(2) * sigma_B;
        b1(i, i+1) = 1i / sqrt(2) * sigma_B;
        b1(i+1, i) = 1 / sqrt(2) * sigma_B;
    end
end
Sigma_u = b1;
% numerical integral
rd = zeros(Dim_U,N);
rd(1:end-2, :) = randn(Dim_Ug * 2 + Dim_UB, N); 
a0 = zeros(Dim_U, 1);
for i = 2:N
    t = i*dt;
    a0(1:2:end-3) = 0;
    a0(2:2:end-2) = 0;
    a0(end-1) = 0;
    a0(end) = 0; 
    u_hat(:,i) = u_hat(:,i-1) + (L_u * u_hat(:,i-1) + a0) * dt + Sigma_u * sqrt(dt) * rd(:, i);
end

% reconstruction
Dim_Grid = 50;
[xx,yy] = meshgrid(linspace(-pi,pi,Dim_Grid), linspace(-pi,pi,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
figure(3) 
for i = 2:2:round(T/dt/100)
    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    xlim([-pi, pi])
    ylim([-pi, pi])
    box on    
    title(['t = ', num2str(dt*100*(i-1))])
    pause(0.1);
    
end


figure
for i = 1:4
    subplot(2,2,i)
    if i == 1
        plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
        title(['(a) gravity mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 2
        plot(dt:dt:N*dt, u_hat(Dim_Ug*2+1,:), 'b', 'linewidth',2)
        title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 3
        plot(dt:dt:N*dt, u_hat(end-1,:), 'b', 'linewidth',2)
        title('(c) Zonal background flow','fontsize',14)
    elseif i == 4
        plot(dt:dt:N*dt, u_hat(end,:), 'b', 'linewidth',2)
        title('(d) Meridional background flow','fontsize',14)
    end
    set(gca,'fontsize',12)
    box on
    xlabel('t')
end