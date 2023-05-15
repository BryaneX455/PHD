% Floe model data assimilation

% dimension of the underlying flow field
Dim_U = length(u_hat(:,1));

% dimension of the unobserved variables v_c, v_o and u_hat
Dim_Y = 3*L + Dim_U;

% dimension of the observed variables [x; y]
Dim_X = 2*L;

% define constant matrices in data assimilation
invBoB = 1 / sigma_x / sigma_x * eye(2*L); % inverse of the square of the observational noise

b1 = [zeros(3*L, 3*L),zeros(3*L, Dim_U);% the noise matrix in the unobserved processes
    zeros(Dim_U, 3*L), Sigma_u];

A1 = eye(2*L, 3*L + Dim_U); % matrix A1 in the observed processes
A0 = zeros(2*L,1); % vector A0 in the observed processes
thickness_min = min(ones(L,1) * thickness, thickness' * ones(1,L));
% save the posterior mean and posterior covariance
gamma_mean_trace = zeros(Dim_Y,N); % posterior mean vector
gamma_cov_trace = zeros(Dim_Y,N); % only same the **diagonal entries** in the posterior covariance matrix
gamma_cov_save = zeros(6,N);
% initial value (can be arbitary)
gamma_mean0 = zeros(Dim_Y,1);
gamma_cov0 = eye(Dim_Y)*0.01; 

gamma_mean_trace(:,1) = gamma_mean0;
gamma_cov_trace(:,1) = diag(gamma_cov0);


% data assimilation
for i = 2:N
    if mod(i,1000) == 0
        disp(i*dt)
    end
    % observational operator 
    x_loc = [x(:,i-1),y(:,i-1)];
    
    G1 = (beta_l./ I  * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) .* (ones(L,1) * (1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:))))/2; % Fourier bases for ocean induced rotation
    G2 = (alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(1,:))); % Fourier bases for u
    G3 = (alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    % computing the difference between the locations in the Lagrangian
    % tracers; need to consider the cases near the boundaries 
    diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + 2*pi; diff_x3 = x(:,i) - x(:,i-1) - 2*pi;  
    diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + 2*pi; diff_y3 = y(:,i) - y(:,i-1) - 2*pi;  
    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    diff_xy = [diff_x; diff_y];
    
    % matrix a0
    % Fmu part
    u_ocn_xf = gamma_mean0(109:184,1).* transpose(rk(1,:));
    u_ocn_yf = gamma_mean0(109:184,1).* transpose(rk(2,:));
    vxf = gamma_mean0(1:36,1);
    vyf = gamma_mean0(37:72,1);
    omgf = gamma_mean0(73:108,1);
    Fm2 = alpha_l ./ m .* (exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf) .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf);
    Fm3 = alpha_l ./ m .* (exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf) .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf);
    Fm1 = beta_l ./ I .* (exp(1i * x_loc * kk) * ( gamma_mean0(109:184,1) .* transpose( 1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:) ) )/2 - omgf) .* abs(exp(1i * x_loc * kk) * ( gamma_mean0(109:184,1) .* transpose( 1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:) ) )/2 - omgf);
    Fm1 = zeros(36,1);
    % Jacobian Nine Block Matrix
    B11 = - 2 * alpha_l ./ m .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf) ;
    B22 = - 2 * alpha_l ./ m .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf) ;
    B33 = - 2 * beta_l ./ I .* abs( exp(1i * x_loc * kk) * ( gamma_mean0(109:184,1) .* transpose( 1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:) ) )/2 - omgf );
    B = diag([B11;B22;B33]);
    B14 = -B11 .* exp(1i * x_loc * kk);
    B15 = -B22 .* exp(1i * x_loc * kk);
    B16 = -B33 .* exp(1i * x_loc * kk);
    BR = [B14;B15;B16];
    Jac = [B,BR;
           zeros(Dim_U, L), zeros(Dim_U, L), zeros(Dim_U, L), L_u];
    
    
    F_u = zeros(Dim_U, 1);
    t = i*dt;
    F_u(1:2:end-3) = 0; 
    F_u(2:2:end-2) = 0; 
    F_u(end-1) = 0;
    F_u(end) = 0; 
    
%     a0 = [Fm2 - Jac(1:36, :) * gamma_mean0; 
%           Fm3 - Jac(37:72, :) * gamma_mean0; 
%           Fm1; 
%           F_u];
%     

    a0 = [Fm2; 
          Fm3; 
          Fm1; 
          F_u];
      
    % matrix a1
%     a1 = [zeros(L,L)*1, zeros(L, L)*1, zeros(L, L)*1, G2;
%           zeros(L, L)*1, zeros(L, L)*1, zeros(L, L)*1, G3;
%           zeros(L, L)*1, zeros(L, L)*1, zeros(L, L)*1, G1;
%          zeros(Dim_U, L), zeros(Dim_U, L), zeros(Dim_U, L), L_u];
    a1 = Jac;
    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A0 * dt - A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     

    % save the posterior statistics
    gamma_mean_trace(:,i) = gamma_mean;
    gamma_cov_trace(:,i) = diag(gamma_cov);

    % update
    gamma_mean0 = gamma_mean;
    gamma_cov0 = gamma_cov;
end

RMSE = zeros(1,Dim_U-2);
PC = zeros(1,Dim_U-2);
for i = 1:Dim_U-2
    RMSE(i) = sqrt(sum((real(u_hat(i,:))-real(gamma_mean_trace(3*L+i,:))).^2)/N);
    corrtemp = corrcoef(real(u_hat(i,:)),real(gamma_mean_trace(3*L+i,:)));
    PC(i) = corrtemp(1,2);
end
RMSE2 = zeros(1,3*L);
PC2 = zeros(1,3*L);
for i = 1:L
    RMSE2(i) = sqrt(sum((real(v_total_x(i,:))-real(gamma_mean_trace(i,:))).^2)/N);
    corrtemp = corrcoef(real(v_total_x(i,:)),real(gamma_mean_trace(i,:)));
    PC2(i) = corrtemp(1,2);
end
for i = 1:L
    RMSE2(L+i) = sqrt(sum((real(v_total_y(i,:))-real(gamma_mean_trace(L+i,:))).^2)/N);
    corrtemp = corrcoef(real(v_total_y(i,:)),real(gamma_mean_trace(L+i,:)));
    PC2(L+i) = corrtemp(1,2);
end
for i = 1:L
    RMSE2(2*L+i) = sqrt(sum((real(omega(i,:))-real(gamma_mean_trace(2*L+i,:))).^2)/N);
    corrtemp = corrcoef(real(omega(i,:)),real(gamma_mean_trace(2*L+i,:)));
    PC2(2*L+i) = corrtemp(1,2);
end

figure
subplot(2,2,1)
plot(RMSE,'b')
box on
set(gca,'fontsize',12)
title('RMSE')
subplot(2,2,2)
plot(PC,'b')
box on
set(gca,'fontsize',12)
title('PC')
subplot(2,2,3)
plot(RMSE2,'b')
box on
set(gca,'fontsize',12)
title('RMSE2')
subplot(2,2,4)
plot(PC2,'b')
box on
set(gca,'fontsize',12)
title('PC2')



% plot the true and the recovered velocity fields
figure
for i = 1:4
    subplot(2,2,i)
    if i == 1
        hold on
        plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, gamma_mean_trace(3*L+1,:), 'r', 'linewidth',2)
        title(['(a) gravity mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 2
        hold on
        plot(dt:dt:N*dt, u_hat(Dim_Ug*2+1,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, gamma_mean_trace(3*L+Dim_Ug*2+1,:), 'r', 'linewidth',2)
        title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
    elseif i == 3
        hold on
        plot(dt:dt:N*dt, u_hat(end-1,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, gamma_mean_trace(end-1,:), 'r', 'linewidth',2)
        title('(c) Zonal background flow','fontsize',14)
    elseif i == 4
        hold on
        plot(dt:dt:N*dt, u_hat(end,:), 'b', 'linewidth',2)
        plot(dt:dt:N*dt, gamma_mean_trace(end,:), 'r', 'linewidth',2)
        title('(d) Meridional background flow','fontsize',14)
    end
    set(gca,'fontsize',12)
    box on
    xlabel('t')
end

figure
for i = 1:3
    for j = 1:2
        subplot(3,2,(i-1)*2+j)
        if i == 1
            hold on
            plot(dt:dt:N*dt, v_total_x(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' trans velocity in x'],'fontsize',14)
        elseif i == 2
            hold on
            plot(dt:dt:N*dt, v_total_y(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(L+j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' trans velocity in y'],'fontsize',14)
        elseif i == 3
            hold on
            plot(dt:dt:N*dt, omega(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(2*L+j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' angular velocity'],'fontsize',14)    
        end
        set(gca,'fontsize',12)
        box on
        xlabel('t')
    end
end
figure
for i = 1:3
    for j = 7:8
        subplot(3,2,(i-1)*2+j-6)
        if i == 1
            hold on
            plot(dt:dt:N*dt, v_total_x(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' translational velocity in x'],'fontsize',14)
        elseif i == 2
            hold on
            plot(dt:dt:N*dt, v_total_y(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(L+j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' translational velocity in y'],'fontsize',14)
        elseif i == 3
            hold on
            plot(dt:dt:N*dt, omega(j,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(2*L+j,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' angular velocity'],'fontsize',14)    
        end
        set(gca,'fontsize',12)
        box on
        xlabel('t')
    end
end