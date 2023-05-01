% Floe model data assimilation
%% EM Algorithm Prelim

% initial guess of parameters
Thic_Est = 3.*thickness_in_km;
n1 = length(Thic_Est);
KK = 50;
Param_save = zeros(n1,KK);
Param_save(:,1) = Thic_Est;
Param_truth = thickness_in_km;

%% Data Assimilation Prelim
% dimension of the underlying flow field; unobserved variables v_c, v_o and
% u_hat; observed variables [x; y]
Dim_U = length(u_hat(:,1));Dim_Y = 3*L + Dim_U;Dim_X = 2*L;


% define constant matrices in data assimilation
invBoB = 1 / sigma_x / sigma_x * eye(2*L); % inverse of the square of the observational noise

b1 = [zeros(3*L, 3*L),zeros(3*L, Dim_U);% the noise matrix in the unobserved processes
    zeros(Dim_U, 3*L), Sigma_u];

A1 = eye(2*L, 3*L + Dim_U); % matrix A1 in the observed processes
A0 = zeros(2*L,1); % vector A0 in the observed processes
%% EM Alg
for k = 2:KK
    % nonlinear smoothing 
    gamma_mean_trace = zeros(Dim_Y,N); % posterior mean vector
    gamma_cov_trace = zeros(Dim_Y,N); % only same the **diagonal entries** in the posterior covariance matrix
    gamma_cov_save = zeros(6,N);
    mu_s = zeros(Dim_Y,N); % smoothing mean
    R_s = zeros(Dim_Y,N); % smoothing variance
    C = zeros(Dim_Y,N); % auxiliary matrix in smoothing
    % initial value (can be arbitary)
    gamma_mean0 = zeros(Dim_Y,1);
    gamma_cov0 = eye(Dim_Y)*0.01;
    gamma_mean_trace(:,1) = gamma_mean0;
    gamma_cov_trace(:,1) = diag(gamma_cov0);
    
    % Filtering
    for i = 2:N
        if mod(i,1000) == 0
            disp(i*dt)
        end
        % observational operator 
        x_loc = [x(:,i-1),y(:,i-1)];
        
    
        G1 = (beta_l./ I  * ones(1,Dim_U)) .* (exp(1i * x_loc * kk / 50.0 * 2 * pi) .* (ones(L,1) * (1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:))))/2; % Fourier bases for ocean induced rotation
        G2 = 8.64*(alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk / 50.0 * 2 * pi) .* (ones(L,1) * rk(1,:)));
        G3 = 8.64*(alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk / 50.0 * 2 * pi) .* (ones(L,1) * rk(2,:))); % Fourier bases for v
    
        
        
         % tracers; need to consider the cases near the boundaries 
        diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + 50.0; diff_x3 = x(:,i) - x(:,i-1) - 50.0;  
        diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + 50.0; diff_y3 = y(:,i) - y(:,i-1) - 50.0;  
        diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
        diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
        diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
        diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
        diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
        diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
        diff_xy = [diff_x; diff_y];
     
       
      
        
    
        % matrix a0, all zeros in our case
        F_u = zeros(Dim_U, 1);
        t = i*dt;
        F_u(1:2:end-3) = 0; 
        F_u(2:2:end-2) = 0; 
        F_u(end-1) = 0;
        F_u(end) = 0; 
        
        a0 = [0*ones(36,1); 
              0*ones(36,1); 
              0*ones(36,1); 
              F_u];
        
        % matrix a1
        a1 = [zeros(L,L)*1, zeros(L, L)*1, zeros(L, L)*1, G2;
              zeros(L, L)*1, zeros(L, L)*1, zeros(L, L)*1, G3;
              zeros(L, L)*1, zeros(L, L)*1, zeros(L, L)*1, G1;
             zeros(Dim_U, L), zeros(Dim_U, L), zeros(Dim_U, L), L_u];
        
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


    % Smoothing
    mu_s(:,end) = gamma_mean; % the starting point is the final value of the filtering
    R_s(:,end) = diag(gamma_cov);
    for i = N-1:-1:1
        % matrices and vectors needed in smoothing
        mu_s(:,i) = mu_s(:,i+1) + (- a0 - a1 * mu_s(:,i+1) + b1 * b1' * inv(diag(gamma_cov_trace(:,i+1)))*(gamma_mean_trace(:,i+1) - mu_s(:,i+1)))*dt; % smoothing mean
        R_s(i) = R_s(:,i+1) - dt*( (a1 + b1 * b1' * inv(diag(gamma_cov_trace(:,i+1))))*R_s(:,i+1) + R_s(:,i+1)*(a1' + b1 * b1' * inv(diag(gamma_cov_trace(:,i+1)))) - b1*b1' ); % smoothing variance
    end
    if k == 2
        mu_s_2 = mu_s;
        R_s_2 = R_s;
    end
    if k == 5
        mu_s_5 = mu_s;
        R_s_5 = R_s;
    end
    MRM = zeros(n1,n1);
    MRZ = zeros(n1,1);

    % M step
    for i = 1+10:N-1-10
        X = mu_s(:,i);
        XX = mu_s(:,i)^2 + R_s(i);
        XXj = mu_s(i) * mu_s(i+1) + R_s(i+1) * C(i)';
        
        MRM_temp = [1/sigma_y^2, - X * z_truth(i)/sigma_y^2, 0, 0;
            -X* z_truth(i)/sigma_y^2, XX * z_truth(i)^2/sigma_y^2 + XX * y_truth(i)^2/sigma_z^2, 0, 0;
            0, 0, XX/sigma_x^2, -X/sigma_x^2;
            0, 0, -X /sigma_x^2, 1/sigma_x^2] * dt^2;
        MRZ_temp = [(y_truth(i+1) - y_truth(i) -  (X * y_truth(i) - y_truth(i)) * dt)/ sigma_y^2;
                    - X * (z_truth(i)/sigma_y^2 * (y_truth(i+1) - y_truth(i) + y_truth(i) * dt)) ...
                    + X * (y_truth(i)/sigma_z^2 * (z_truth(i+1) - z_truth(i) + z_truth(i) * dt)) ...
                    + y_truth(i) * z_truth(i)/sigma_y^2 * dt * XX - y_truth(i) * z_truth(i)/sigma_z^2 * dt * XX;
                    -XXj/sigma_x^2 + XX/sigma_x^2 - (y_truth(i)^2 + z_truth(i)^2) /sigma_x^2 * X * dt;
                    (mu_s(i+1) - X + (y_truth(i)^2 + z_truth(i)^2) * dt)/sigma_x^2]*dt; 
                    
        
        MRM = MRM + MRM_temp;
        MRZ = MRZ + MRZ_temp;
    end
    Theta = MRM \ MRZ; % quadratic optimization
end

%% Showing results
% error assignment
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
% plots
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