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
% gamma_cov0 = zeros(Dim_Y);
% gamma_cov0 = cov_temp;
% for i = 1:3*L
%     gamma_cov0(i,i) = 0.005;
% end
% for i = 3*L+1: 3*L+Dim_Ug*2
%     gamma_cov0(i,i) = sigma_g^2/(d_g + sqrt(d_g^2 + L*sigma_x^(-2)*sigma_g^2));
% end
% for i = 3*L+Dim_Ug * 2: 3*L+Dim_Ug * 2 + Dim_UB
%     gamma_cov0(i,i) = sigma_B^2/(d_B + sqrt(d_B^2 + L*sigma_x^(-2)*sigma_B^2));
% end

% data assimilation
for i = 2:N
    if mod(i,1000) == 0
        disp(i*dt)
    end
    % observational operator 
    x_loc = [x(:,i-1),y(:,i-1)];
    
%     G1 = (beta_l./ I  * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) * diag(transpose(1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:))))/2; % Fourier bases for ocean induced rotation
%     G2 = (alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) * diag(transpose(rk(1,:)))); % Fourier bases for u
%     G3 = (alpha_l./ m * ones(1,Dim_U)) .* (exp(1i * x_loc * kk) * diag(transpose(rk(2,:)))); % Fourier bases for v
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
    
    % finding the contact floes
    distance_x1 = abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)'); % computing the distance in the general case inside the domain
    distance_x2 = abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + 2*pi); % dealing with the floes near the boundaries
    distance_x3 = abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' - 2*pi); % dealing with the floes near the boundaries
    distance_x_temp = min(distance_x1, distance_x2); % the distance in the x direction
    distance_x = min(distance_x_temp, distance_x3); % the distance in the x direction    
    distance_y1 = abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)'); % computing the distance in the general case inside the domain
    distance_y2 = abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + 2*pi); % dealing with the floes near the boundaries
    distance_y3 = abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' - 2*pi); % dealing with the floes near the boundaries
    distance_y_temp = min(distance_y1, distance_y2); % the distance in the y direction
    distance_y = min(distance_y_temp, distance_y3); % the distance in the y direction    
    distance = sqrt(distance_x.^2 + distance_y.^2); % the overall distance
    distance_selected = (distance < distance_every_two_floes);
    projection_x = distance_x ./ (distance + eps * eye(L)) .* distance_selected; 
    projection_y = distance_y ./ (distance + eps * eye(L)) .* distance_selected;
    
    equal_temp_x1 = (distance_x == distance_x1) .* distance_selected; % finding the general case inside the domain
    equal_temp_x2 = (distance_x == distance_x2) .* distance_selected; % finding the case near the boundary
    equal_temp_x3 = (distance_x == distance_x3) .* distance_selected; % finding the case near the boundary
    distance_vector_x1_temp = (distance_every_two_floes .* projection_x - distance_x1) .* (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)') ...
        ./ abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + eps * eye(L)) .* equal_temp_x1; % compressive distance in x direction
    % compression length in x direction * unit vector in x (determining the positive or negative sign) * the selected cases
    distance_vector_x2_temp = (distance_every_two_floes .* projection_x - distance_x2) .* (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + 2*pi) ...
        ./ abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + 2*pi).* equal_temp_x2;
    distance_vector_x3_temp = (distance_every_two_floes.* projection_x - distance_x3) .* (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' - 2*pi) ...
        ./ abs(x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' - 2*pi).* equal_temp_x3;
    
    % computing the transverse_area (namely the chord length)
    transverse_area = 1 ./ (distance + eps * eye(L)) .* sqrt(4 * distance.^2 .* radius_max.^2 - (distance.^2 - radius_min.^2 + radius_max.^2).^2) .* distance_selected;
    
    
    normal_force_x = - Elj .* thickness_min .* transverse_area .* (distance_vector_x1_temp + distance_vector_x2_temp + distance_vector_x3_temp); % resistence forcing in the x direction

    
    equal_temp_y1 = (distance_y == distance_y1) .* distance_selected; % finding the general case inside the domain
    equal_temp_y2 = (distance_y == distance_y2) .* distance_selected; % finding the case near the boundary
    equal_temp_y3 = (distance_y == distance_y3) .* distance_selected; % finding the case near the boundary
    distance_vector_y1_temp = (distance_every_two_floes .* projection_y - distance_y1) .* (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)') ...
        ./ abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + eps * eye(L)) .* equal_temp_y1; % compressive distance in y direction
    % compression length in y direction * unit vector in x (determining the positive or negative sign) * the selected cases
    distance_vector_y2_temp = (distance_every_two_floes .* projection_y - distance_y2) .* (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + 2*pi) ...
        ./ abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + 2*pi).* equal_temp_y2;
    distance_vector_y3_temp = (distance_every_two_floes.* projection_y - distance_y3) .* (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' - 2*pi) ...
        ./ abs(y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' - 2*pi).* equal_temp_y3;
    normal_force_y = - Elj .* thickness_min .* transverse_area .* (distance_vector_y1_temp + distance_vector_y2_temp + distance_vector_y3_temp); % resistence forcing in the y direction
    
    % computing the normal and tangential unit vectors
    normal_direction_x = ((x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)') .* equal_temp_x1 +...
        (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + 2*pi) .* equal_temp_x2 +...
        (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' - 2*pi) .* equal_temp_x3) ./ (distance + eps * eye(L));
    normal_direction_y = ((y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)') .* equal_temp_y1 +...
        (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + 2*pi) .* equal_temp_y2 +...
        (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' - 2*pi) .* equal_temp_y3) ./ (distance + eps * eye(L));
    tangential_direction_x = normal_direction_y;
    tangential_direction_y = - normal_direction_x;
%     transverse_area = 1 ./ (distance + eps * eye(L)) .* sqrt(4 * distance.^2 .* radius_max.^2 - (distance.^2 - radius_min.^2 + radius_max.^2).^2) .* distance_selected;
    
    % compute the tangential force due to the contraints
    % the true velocity and angular velocity are both unknown. the filtered
    % mean values are used in the estimation here
    v_diff_x = + real( gamma_mean_trace(1:L,i-1) * ones(1,L) - ones(L,1) * gamma_mean_trace(1:L,i-1)' );
    v_diff_y = + real( gamma_mean_trace(L+1:2*L,i-1) * ones(1,L) - ones(L,1) * gamma_mean_trace(L+1:2*L,i-1)' );

    % computing difference between the angular velocity multiplying the radius at the contact interface
    angular_velocity_diff_x = - (ones(L,1) * radius)' .* normal_direction_y' .* (ones(L,1) * real(gamma_mean_trace(2*L+1:3*L,i-1))')' + (ones(L,1) * radius) .* normal_direction_y .* (ones(L,1) * real(gamma_mean_trace(2*L+1:3*L,i-1))');
    angular_velocity_diff_y =   (ones(L,1) * radius)' .* normal_direction_x' .* (ones(L,1) * real(gamma_mean_trace(2*L+1:3*L,i-1))')' - (ones(L,1) * radius) .* normal_direction_x .* (ones(L,1) * real(gamma_mean_trace(2*L+1:3*L,i-1))');
    
 

%     v_diff_x = + real( v_total_x(:,i-1) * ones(1,L) - ones(L,1) * v_total_x(:,i-1)' );
%     v_diff_y = + real( v_total_y(:,i-1) * ones(1,L) - ones(L,1) * v_total_y(:,i-1)' );
% 
%     % computing difference between the angular velocity multiplying the radius at the contact interface
%     angular_velocity_diff_x = - (ones(L,1) * radius)' .* normal_direction_y' .* (ones(L,1) * real(omega(:,i-1))')' + (ones(L,1) * radius) .* normal_direction_y .* (ones(L,1) * real(omega(:,i-1))');
%     angular_velocity_diff_y =   (ones(L,1) * radius)' .* normal_direction_x' .* (ones(L,1) * real(omega(:,i-1))')' - (ones(L,1) * radius) .* normal_direction_x .* (ones(L,1) * real(omega(:,i-1))');
    
 


    % computing the total tangential force and its projections to the x and
    % y axes
    tangential_force = ((v_diff_x + angular_velocity_diff_x) .* tangential_direction_x + (v_diff_y + angular_velocity_diff_y) .* tangential_direction_y);
    tangential_force_temp = mulj .* sqrt(normal_force_x.^2 + normal_force_y.^2);
    tangential_force_max = min(tangential_force_temp, abs(tangential_force));
    
    
    distance_selected_1 = (tangential_force_max == tangential_force_temp) .* (tangential_force_max ~= 0);
    distance_selected_2 = (tangential_force_max == abs(tangential_force)) .* (tangential_force_max ~= 0);
    
    tangential_force_1 = ( (tangential_force_max == tangential_force_temp) .* tangential_force_temp ) .* sign(tangential_force);
    
    tangential_force_x_1 = Glj .* thickness_min .* transverse_area .* tangential_force_1 .* tangential_direction_x;
    tangential_force_y_1 = Glj .* thickness_min .* transverse_area .* tangential_force_1 .* tangential_direction_y;
    
    rotation_force_a0_part = 1./I .* sum( (ones(L,1) * radius)' .* (normal_direction_x .* tangential_force_y_1 - normal_direction_y .* tangential_force_x_1))';
    % matrix for the tangential force in front of v_o+v_c 
    a1_v1_temp1 = Glj .* thickness_min .* transverse_area .* distance_selected_2 .* (tangential_direction_x.^2);    
    a1_v1_temp2 = diag(sum(-a1_v1_temp1)) + a1_v1_temp1; % the summation is because one floe can be hit by more than one floes; the negative sign is the difference between the two floe velocities
    a1_v2_temp1 = Glj .* thickness_min .* transverse_area .* distance_selected_2 .* (tangential_direction_x .* tangential_direction_y);    
    a1_v2_temp2 = diag(sum(-a1_v2_temp1)) + a1_v2_temp1;    
    a1_v1_temp3 = Glj .* thickness_min .* transverse_area .* distance_selected_2 .* (tangential_direction_x .* tangential_direction_y);    
    a1_v1_temp4 = diag(sum(-a1_v1_temp3)) + a1_v1_temp3; % the summation is because one floe can be hit by more than one floes; the negative sign is the difference between the two floe velocities
    a1_v2_temp3 = Glj .* thickness_min .* transverse_area .* distance_selected_2 .* (tangential_direction_y.^2);    
    a1_v2_temp4 = diag(sum(-a1_v2_temp3)) + a1_v2_temp3;   
    
    a1_omega_temp1 = - Glj .* transverse_area .* thickness_min .* distance_selected_2 .* ( ones(L,1) * radius ) .* ( normal_direction_y .* (tangential_direction_x.^2) - normal_direction_x .* (tangential_direction_x .* tangential_direction_y) ); 
    a1_omega_temp2 = diag(sum(-a1_omega_temp1)) + a1_omega_temp1;   
    a1_omega_temp3 = - Glj .* transverse_area .* thickness_min .* distance_selected_2 .* ( ones(L,1) * radius ) .* ( normal_direction_y .* (tangential_direction_x .* tangential_direction_y) - normal_direction_x .* (tangential_direction_y.^2) ); 
    a1_omega_temp4 = diag(sum(-a1_omega_temp3)) + a1_omega_temp3;  
    
    a1_11 = (1./m * ones(1,L)) .* a1_v1_temp2  - diag(alpha_l./ m);
    a1_12 = (1./m * ones(1,L)) .* a1_v2_temp2 ;
    a1_13 = (1./m * ones(1,L)) .* a1_omega_temp2 ;
%     if mod(i,100) == 0
%     (1./m * ones(1,L)) .* a1_v1_temp2
%     pause
%     end
    a1_21 = (1./m * ones(1,L)) .* a1_v1_temp4;
    a1_22 = (1./m * ones(1,L)) .* a1_v2_temp4  - diag(alpha_l./ m);    
    a1_23 = (1./m * ones(1,L)) .* a1_omega_temp4 ;
    
    a1_31_temp  = - normal_direction_y' .* a1_v1_temp1 + normal_direction_x' .* a1_v1_temp3;
    a1_31_temp2 = diag(sum(- normal_direction_y .* (-a1_v1_temp1) + normal_direction_x .* (-a1_v1_temp3)));
    a1_31 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* (a1_31_temp + a1_31_temp2);
    
    a1_32_temp  = - normal_direction_y' .* a1_v2_temp1 + normal_direction_x' .* a1_v2_temp3;
    a1_32_temp2 = diag(sum(- normal_direction_y .* (-a1_v2_temp1) + normal_direction_x .* (-a1_v2_temp3)));
    a1_32 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* (a1_32_temp + a1_32_temp2);
    
    a1_33_temp  = - normal_direction_y' .* a1_omega_temp1 + normal_direction_x' .* a1_omega_temp3;
    a1_33_temp2 = diag(sum(- normal_direction_y .* (-a1_omega_temp1) + normal_direction_x .* (-a1_omega_temp3)));
    a1_33 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* (a1_33_temp + a1_33_temp2) - diag(beta_l./I);

%     a1_31 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* ( - normal_direction_y .* a1_v1_temp2 + normal_direction_x .* a1_v1_temp4 ) ;
%     a1_32 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* ( - normal_direction_y .* a1_v2_temp2 + normal_direction_x .* a1_v2_temp4 ) ;
%     a1_33 = (1./I * ones(1,L)) .* (radius' * ones(1,L)) .* ( - normal_direction_y .* a1_omega_temp2 + normal_direction_x .* a1_omega_temp4 ) - diag(beta_l./I);
    

    % matrix a0
    F_u = zeros(Dim_U, 1);
    t = i*dt;
    F_u(1:2:end-3) =0; 0.4+0.4*1i;%f_amp * exp(1i * f_phase * t) * ones(Dim_Ug + Dim_UB/2, 1);
    F_u(2:2:end-2) =0; 0.4-0.4*1i;%f_amp * exp(- 1i * f_phase * t) * ones(Dim_Ug + Dim_UB/2, 1);
    F_u(end-1) = 0;f_amp * cos(f_phase * t) + f_x_b;0;
    F_u(end) = 0;f_amp * sin(f_phase * t) + f_y_b; 0;   
    
    a0 = [1./m .* sum(normal_force_x)' + 1./m .* sum(tangential_force_x_1)'; 
          1./m .* sum(normal_force_y)' + 1./m .* sum(tangential_force_y_1)'; 
          rotation_force_a0_part; 
          F_u];
    
    % matrix a1
    a1 = [a1_11, a1_12, a1_13, G2;
          a1_21, a1_22, a1_23, G3;
          a1_31, a1_32, a1_33, G1;
         zeros(Dim_U, L), zeros(Dim_U, L), zeros(Dim_U, L), L_u];
    
    % run the data assimilation for posterior mean and posterior covariance
    gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (diff_xy - A0*dt-A1 * gamma_mean0 * dt);
    gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
%     gamma_mean(1:3*L) = real(gamma_mean(1:3*L));
% gamma_mean(1:12) 
% pause
% gamma_cov_save(1,i) = gamma_cov(1,241);
% gamma_cov_save(2,i) = gamma_cov(1,259);
% gamma_cov_save(3,i) = gamma_cov(2,241);
% gamma_cov_save(4,i) = gamma_cov(2,259);
% gamma_cov_save(5,i) = gamma_cov(2*L+1,241);
% gamma_cov_save(6,i) = gamma_cov(2*L+1,259);
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