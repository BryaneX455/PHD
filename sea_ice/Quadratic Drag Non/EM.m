% Floe model data assimilationc
% dimension of the underlying flow field
obs_num = 3; emit = 2;
Dim_U = length(u_hat(:,1));

% dimension of the unobserved variables v_c, v_o and u_hat
Dim_Y = 3*L + Dim_U;

% dimension of the observed variables [x; y]
Dim_X = obs_num*L;

% define constant matrices in data assimilation
invBoB = 1 / sigma_x / sigma_x * eye(obs_num*L); % inverse of the square of the observational noise

b1 = [diag(ones(L,1)*sigma_x),zeros(L,2*L+Dim_U);
    zeros(L,L),diag(ones(L,1)*sigma_x),zeros(L,L+Dim_U);
    zeros(L,2*L),diag(ones(L,1)*sigma_x),zeros(L,Dim_U);
    zeros(Dim_U, 3*L), Sigma_u];

A1 = eye(obs_num*L, 3*L + Dim_U); % matrix A1 in the observed processes
A0 = zeros(obs_num*L,1); % vector A0 in the observed processes
thickness_min = min(ones(L,1) * thickness, thickness' * ones(1,L));
% Initial Guess
theta = ones(L,1)/2;
for k = 1:emit-1
    % save the posterior mean and posterior covariance
    gamma_mean_trace = zeros(Dim_Y,N); % posterior mean vector
    gamma_cov_trace = zeros(Dim_Y,N); % only same the **diagonal entries** in the posterior covariance matrix
    gamma_cov_save = zeros(6,N);
    % initial value (can be arbitary)
    gamma_mean0 = zeros(Dim_Y,1);
    gamma_cov0 = eye(Dim_Y)*0.01; 

    gamma_mean_trace(:,1) = gamma_mean0;
    gamma_cov_trace(:,1) = diag(gamma_cov0);


    Thickness = zeros(L,emit);
    % data assimilation
    disp(theta);
    for i = 2:N
        if mod(i,1000) == 0
            disp(i*dt)
        end
        % observational operator 
        x_loc = [x(:,i-1),y(:,i-1)];

        % computing the difference between the locations in the Lagrangian tracers; need to consider the cases near the boundaries 
        diff_x1 = x(:,i) - x(:,i-1); diff_x2 = x(:,i) - x(:,i-1) + 2*pi; diff_x3 = x(:,i) - x(:,i-1) - 2*pi;  
        diff_y1 = y(:,i) - y(:,i-1); diff_y2 = y(:,i) - y(:,i-1) + 2*pi; diff_y3 = y(:,i) - y(:,i-1) - 2*pi;  
        diff_omg1 = Omg(:,i) - Omg(:,i-1); diff_omg2 = Omg(:,i) - Omg(:,i-1) + 2*pi; diff_omg3 = Omg(:,i) - Omg(:,i-1) - 2*pi; % new added angular velocity part

        diff_xtemp = min(abs(diff_x1), abs(diff_x2)); diff_x_index = min(abs(diff_x3), diff_xtemp);
        diff_ytemp = min(abs(diff_y1), abs(diff_y2)); diff_y_index = min(abs(diff_y3), diff_ytemp);
        diff_omgtemp = min(abs(diff_omg1), abs(diff_omg2)); diff_omg_index = min(abs(diff_omg3), diff_omgtemp); % added angular part

        diff_x1_index = (diff_x_index == abs(diff_x1)); diff_x2_index = (diff_x_index == abs(diff_x2)); diff_x3_index = (diff_x_index == abs(diff_x3)); 
        diff_y1_index = (diff_y_index == abs(diff_y1)); diff_y2_index = (diff_y_index == abs(diff_y2)); diff_y3_index = (diff_y_index == abs(diff_y3)); 
        diff_omg1_index = (diff_omg_index == abs(diff_omg1));diff_omg2_index = (diff_omg_index == abs(diff_omg2));diff_omg3_index = (diff_omg_index == abs(diff_omg3));

        diff_x = diff_x1 .* diff_x1_index + diff_x2 .* diff_x2_index + diff_x3 .* diff_x3_index;
        diff_y = diff_y1 .* diff_y1_index + diff_y2 .* diff_y2_index + diff_y3 .* diff_y3_index;
        diff_omg = diff_omg1 .* diff_omg1_index + diff_omg2 .* diff_omg2_index + diff_omg3 .* diff_omg3_index;
        diff_xy = [diff_x; diff_y; diff_omg];
        % diff_xy = [diff_x; diff_y];

        % matrix a0
        % Fmu part
        u_ocn_xf = gamma_mean0(3*L+1:end,1).* transpose(rk(1,:));
        u_ocn_yf = gamma_mean0(3*L+1:end,1).* transpose(rk(2,:));
        vxf = gamma_mean0(1:L,1);
        vyf = gamma_mean0(L+1:2*L,1);
        omgf = gamma_mean0(2*L+1:3*L,1);
        Fm2 = theta .* (exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf) .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf);
        Fm3 = theta .* (exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf) .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf);
        Fm1 = theta .* (exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf) .* abs(exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf);

        % Jacobian 
        JR11 = - 2 * (theta) .* abs( exp(1i * x_loc * kk) * (u_ocn_xf)  - vxf);
        JR22 = - 2 * (theta) .* abs( exp(1i * x_loc * kk) * (u_ocn_yf)  - vyf);
        JR33 = - 2 * (theta) .* abs( exp(1i * x_loc * kk) * ( gamma_mean0(3*L+1:end,1) .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgf );
        JR_Vel = diag([JR11;JR22;JR33]);

        % V1
        JR14 = - JR11 .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(1,:))); 
        JR24 = - JR22 .* (exp(1i * x_loc * kk) .* (ones(L,1) * rk(2,:)));
        JR34 = - JR33 .* (exp(1i * x_loc * kk) .* (ones(L,1) * ( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:)))/2);


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
        gamma_mean = gamma_mean0 + Fm * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A1 * gamma_mean0 * dt);
        gamma_cov = gamma_cov0 + (Jac * gamma_cov0 + gamma_cov0 * Jac' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
        % save the posterior statistics
        gamma_mean_trace(:,i) = gamma_mean;
        gamma_cov_trace(:,i) = diag(gamma_cov);
        % update
        gamma_mean0 = gamma_mean;
        gamma_cov0 = gamma_cov;

    end


    % MRM
    MRM = zeros(L,L);
    MRZ = zeros(L,1);
    Rinv = [diag(1./(ones(108,1) * sigma_x).^2),zeros(108,184);zeros(108,108),diag(1./(ones(108,1) * sigma_x).^2),zeros(108,76);zeros(74,216),diag(1./diag(Sigma_u(1:74,1:74)).^2),zeros(74,2);zeros(2,292)]; % With Sigma_u
    for j = 2:N-1
        ujp1 = gamma_mean_trace(1:L,j+1);
        uj = gamma_mean_trace(1:L,j);
        vjp1 = gamma_mean_trace(L+1:2*L,j+1);
        vj = gamma_mean_trace(L+1:2*L,j);
        omgjp1 = gamma_mean_trace(2*L+1:3*L,j+1);
        omgj = gamma_mean_trace(2*L+1:3*L,j);
        uvjp1 = gamma_mean_trace(3*L+1:end,j+1);
        uvj = gamma_mean_trace(3*L+1:end,j);
        Zjp1 = [x(:,j+1);y(:,j+1);Omg(:,j+1);ujp1;vjp1;omgjp1;uvjp1];
        Zj = [x(:,j);y(:,j);Omg(:,j);uj;vj;omgj;uvj];
        Cj = [uj;vj;omgj;zeros(108,1);F_u + L_u * uvj];
        FM2j = diag(theta .* (exp(1i * x_loc * kk) * (uvj.* transpose(rk(1,:)))  - uj) .* abs( exp(1i * x_loc * kk) * (uvj.* transpose(rk(1,:)))  - uj));
        FM3j = diag(theta .* (exp(1i * x_loc * kk) * (uvj.* transpose(rk(2,:)))  - vj) .* abs( exp(1i * x_loc * kk) * (uvj.* transpose(rk(2,:)))  - vj));
        FM1j = diag(theta .* (exp(1i * x_loc * kk) * ( uvj .* transpose(1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgj) .* abs(exp(1i * x_loc * kk) * ( uvj .* transpose( 1i * rk(2,:) .* kk(1,:) - 1i * rk(1,:) .* kk(2,:) ) )/2 - omgj));
        Mj = [zeros(108,36);FM2j;FM3j;FM1j;zeros(76,36)];
        MRMTemp = Mj' * Rinv * Mj * dt;
        MRZTemp = Mj' * Rinv * (Zjp1 - Zj - Cj*dt);
        MRM = MRM + MRMTemp;
        MRZ = MRZ + MRZTemp;
    end
    theta = MRM\MRZ;
    Thick_Est = d_o./theta;
    Thickness(:,k+1) = Thick_Est;
    disp(theta);
end


% % plot
% 
% RMSE = zeros(1,Dim_U-2);
% PC = zeros(1,Dim_U-2);
% for i = 1:Dim_U-2
%     RMSE(i) = sqrt(sum((real(u_hat(i,:))-real(gamma_mean_trace(3*L+i,:))).^2)/N);
%     corrtemp = corrcoef(real(u_hat(i,:)),real(gamma_mean_trace(3*L+i,:)));
%     PC(i) = corrtemp(1,2);
% end
% RMSE2 = zeros(1,3*L);
% PC2 = zeros(1,3*L);
% for i = 1:L
%     RMSE2(i) = sqrt(sum((real(v_total_x(i,:))-real(gamma_mean_trace(i,:))).^2)/N);
%     corrtemp = corrcoef(real(v_total_x(i,:)),real(gamma_mean_trace(i,:)));
%     PC2(i) = corrtemp(1,2);
% end
% for i = 1:L
%     RMSE2(L+i) = sqrt(sum((real(v_total_y(i,:))-real(gamma_mean_trace(L+i,:))).^2)/N);
%     corrtemp = corrcoef(real(v_total_y(i,:)),real(gamma_mean_trace(L+i,:)));
%     PC2(L+i) = corrtemp(1,2);
% end
% for i = 1:L
%     RMSE2(2*L+i) = sqrt(sum((real(omega(i,:))-real(gamma_mean_trace(2*L+i,:))).^2)/N);
%     corrtemp = corrcoef(real(omega(i,:)),real(gamma_mean_trace(2*L+i,:)));
%     PC2(2*L+i) = corrtemp(1,2);
% end
% 
% figure
% subplot(2,2,1)
% plot(RMSE,'b')
% box on
% set(gca,'fontsize',12)
% title('RMSE')
% subplot(2,2,2)
% plot(PC,'b')
% box on
% set(gca,'fontsize',12)
% title('PC')
% subplot(2,2,3)
% plot(RMSE2,'b')
% box on
% set(gca,'fontsize',12)
% title('RMSE2')
% subplot(2,2,4)
% plot(PC2,'b')
% box on
% set(gca,'fontsize',12)
% title('PC2')
% 
% 
% 
% % plot the true and the recovered velocity fields
% figure
% for i = 1:4
%     subplot(2,2,i)
%     if i == 1
%         hold on
%         plot(dt:dt:N*dt, u_hat(1,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, gamma_mean_trace(3*L+1,:), 'r', 'linewidth',2)
%         title(['(a) gravity mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
%     elseif i == 2
%         hold on
%         plot(dt:dt:N*dt, u_hat(Dim_Ug*2+1,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, gamma_mean_trace(3*L+Dim_Ug*2+1,:), 'r', 'linewidth',2)
%         title(['(b) GB mode ( ', num2str(kk(1,1)),' , ', num2str(kk(2,1)), ' )'],'fontsize',14)
%     elseif i == 3
%         hold on
%         plot(dt:dt:N*dt, u_hat(end-1,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, gamma_mean_trace(end-1,:), 'r', 'linewidth',2)
%         title('(c) Zonal background flow','fontsize',14)
%     elseif i == 4
%         hold on
%         plot(dt:dt:N*dt, u_hat(end,:), 'b', 'linewidth',2)
%         plot(dt:dt:N*dt, gamma_mean_trace(end,:), 'r', 'linewidth',2)
%         title('(d) Meridional background flow','fontsize',14)
%     end
%     set(gca,'fontsize',12)
%     box on
%     xlabel('t')
% end
% 
% figure
% for i = 1:3
%     for j = 1:2
%         subplot(3,2,(i-1)*2+j)
%         if i == 1
%             hold on
%             plot(dt:dt:N*dt, v_total_x(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' trans velocity in x'],'fontsize',14)
%         elseif i == 2
%             hold on
%             plot(dt:dt:N*dt, v_total_y(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(L+j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' trans velocity in y'],'fontsize',14)
%         elseif i == 3
%             hold on
%             plot(dt:dt:N*dt, omega(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(2*L+j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' angular velocity'],'fontsize',14)    
%         end
%         set(gca,'fontsize',12)
%         box on
%         xlabel('t')
%     end
% end
% figure
% for i = 1:3
%     for j = 7:8
%         subplot(3,2,(i-1)*2+j-6)
%         if i == 1
%             hold on
%             plot(dt:dt:N*dt, v_total_x(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' translational velocity in x'],'fontsize',14)
%         elseif i == 2
%             hold on
%             plot(dt:dt:N*dt, v_total_y(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(L+j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' translational velocity in y'],'fontsize',14)
%         elseif i == 3
%             hold on
%             plot(dt:dt:N*dt, omega(j,:), 'b', 'linewidth',2)
%             plot(dt:dt:N*dt, gamma_mean_trace(2*L+j,:), 'r', 'linewidth',2)
%             title(['Floe # ', num2str(j),' angular velocity'],'fontsize',14)    
%         end
%         set(gca,'fontsize',12)
%         box on
%         xlabel('t')
%     end
% end