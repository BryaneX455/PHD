% Floe model data assimilationc
% dimension of the underlying flow field
obs_num = 2;
FL = 36;
Dim_U = length(u_hat(:,1));

% dimension of the unobserved variables v_c, v_o and u_hat
Dim_Y = 3*L + Dim_U;

% dimension of the observed variables [x; y]
Dim_X = obs_num*FL;

% define constant matrices in data assimilation
invBoB = 1 / sigma_x / sigma_x * eye(obs_num*FL); % inverse of the square of the observational noise

b1 = [diag(ones(L,1)*sigma_x),zeros(L,2*L+Dim_U);
    zeros(L,L),diag(ones(L,1)*sigma_x),zeros(L,L+Dim_U);
    zeros(L,2*L),diag(ones(L,1)*sigma_x),zeros(L,Dim_U);
    zeros(Dim_U, 3*L), Sigma_u];

A1 = eye(obs_num*FL, 3*L + Dim_U); % matrix A1 in the observed processes
% A0 = zeros(obs_num*FL,1); % vector A0 in the observed processes
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

% x(FL+1:end,:) = 0;
% y(FL+1:end,:) = 0;
% Omg(FL+1:end,:) = 0;
% data assimilation
for i = 2:N
    if mod(i,1000) == 0
        disp(i*dt)
    end
    % reduce floe number in the vector
    % observational operator 
    x_loc = [x(1:L,i-1),y(1:L,i-1)];
    
    % computing the difference between the locations in the Lagrangian
    % tracers; need to consider the cases near the boundaries 
    diff_x1 = x(1:L,i) - x(1:L,i-1); diff_x2 = x(1:L,i) - x(1:L,i-1) + 2*pi; diff_x3 = x(1:L,i) - x(1:L,i-1) - 2*pi;  
    diff_y1 = y(1:L,i) - y(1:L,i-1); diff_y2 = y(1:L,i) - y(1:L,i-1) + 2*pi; diff_y3 = y(1:L,i) - y(1:L,i-1) - 2*pi;  
    diff_omg1 = Omg(1:L,i) - Omg(1:L,i-1); diff_omg2 = Omg(1:L,i) - Omg(1:L,i-1) + 2*pi; diff_omg3 = Omg(1:L,i) - Omg(1:L,i-1) - 2*pi; % new added angular velocity part

    diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
    diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
    diff_omgtemp = min(abs(diff_omg1), abs(diff_omg2)); diff_omg_index = min(abs(diff_omg3), diff_omgtemp); % added angular part

    diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
    diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
    diff_omg1_index = (diff_omg_index == abs(diff_omg1));diff_omg2_index = (diff_omg_index == abs(diff_omg2));diff_omg3_index = (diff_omg_index == abs(diff_omg3));

    diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
    diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
    diff_omg = diff_omg1 .* diff_omg1_index + diff_omg2 .* diff_omg2_index + diff_omg3 .* diff_omg3_index;
%     diff_xy = [diff_x(1:FL); diff_y(1:FL); diff_omg(1:FL)];
    diff_xy = [diff_x; diff_y];
    
    % matrix a0
    % Fmu part
    u_ocn_xf = gamma_mean0(3*L+1:end,1).* transpose(rk(1,:));
    u_ocn_yf = gamma_mean0(3*L+1:end,1).* transpose(rk(2,:));
    vxf = gamma_mean0(1:L,1);
    vyf = gamma_mean0(L+1:2*L,1);
    omgf = gamma_mean0(2*L+1:3*L,1);
    Fm2 = alpha_l(1:L,1) ./ m(1:L,1) .* (exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf) .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf) + sqrt(dt) * sigma_x * randn(L,1);
    Fm3 = alpha_l(1:L,1) ./ m (1:L,1).* (exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf) .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf) + sqrt(dt) * sigma_x * randn(L,1);
    Fm1 = beta_l(1:L,1) ./ I(1:L,1) .* (exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf) .* abs(exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf) + sqrt(dt) * sigma_x * randn(L,1);

    % Jacobian 
    JR11 = - 2 * (alpha_l(1:L,1) ./ m(1:L,1)) .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf);
    JR22 = - 2 * (alpha_l(1:L,1) ./ m(1:L,1)) .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf);
    JR33 = - 2 * (beta_l(1:L,1) ./ I(1:L,1)) .* abs( exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf );
    JR_Vel = diag([JR11;JR22;JR33]);

    % V1
    JR14 = - JR11 .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(1,:))); 
    JR24 = - JR22 .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(2,:)));
    JR34 = - JR33 .* (exp(1i * x_loc * kk) .* (ones(L,1) * (1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:)))/2);
   

    JR_Ocn = [JR14;JR24;JR34];
    Jac = [JR_Vel,JR_Ocn;
           zeros(Dim_U, L), zeros(Dim_U, L), zeros(Dim_U, L), L_u];
    
    
    F_u = zeros(Dim_U, 1);
    t = i*dt;
    F_u(1:2:end-3) = 0; 
    F_u(2:2:end-2) = 0; 
    F_u(end-1) = 0;
    F_u(end) = 0; 
    
    a0 = [Fm2 - Jac(1:L, :) * gamma_mean0; 
          Fm3 - Jac(L+1:2*L, :) * gamma_mean0; 
          Fm1 - Jac(2*L+1:3*L, :) * gamma_mean0; 
          F_u];
      
    Fm = [Fm2; 
          Fm3; 
          Fm1; 
          F_u + L_u * gamma_mean0(3*L+1:end,1)];

      

    a1 = Jac;
    % run the data assimilation for posterior mean and posterior covariance
    % gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A1 * gamma_mean0 * dt);
    % gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
    % The two gamma_mean(s) are the same, refer to the analytical part
    gamma_mean = gamma_mean0 + Fm * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (Jac * gamma_cov0 + gamma_cov0 * Jac' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
    % save the posterior statistics
    gamma_mean_trace(:,i) = gamma_mean;
%     gamma_cov_trace(:,i) = diag(gamma_cov);
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
            plot(dt:dt:N*dt, v_total_x(25,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(25,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' translational velocity in x'],'fontsize',14)
        elseif i == 2
            hold on
            plot(dt:dt:N*dt, v_total_y(25,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(L+25,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' translational velocity in y'],'fontsize',14)
        elseif i == 3
            hold on
            plot(dt:dt:N*dt, omega(25,:), 'b', 'linewidth',2)
            plot(dt:dt:N*dt, gamma_mean_trace(2*L+25,:), 'r', 'linewidth',2)
            title(['Floe # ', num2str(j),' angular velocity'],'fontsize',14)    
        end
        set(gca,'fontsize',12)
        box on
        xlabel('t')
    end
end