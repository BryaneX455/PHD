% Floe model main file
rng(11);
T = 10; % total time length
dt = 0.001; % time step
N = T/dt; % total time points in time
x = zeros(L,N);
y = zeros(L,N);
x(:,1) = Location(1,:)'; % Lagrangian tracers (floe center location x) evolving with time
y(:,1) = Location(2,:)'; % Lagrangian tracers (floe center location y) evolving with time
vc_x = zeros(L,N); % velocity component x corresponding to the contact forces
vc_y = zeros(L,N); % velocity component y corresponding to the contact forces
vo_x = zeros(L,N); % velocity component x corresponding to the ocean forces
vo_y = zeros(L,N); % velocity component y corresponding to the ocean forces
sigma_x = 0.1/2; % random noise in the tracer equation
omega = zeros(L,N); % angular velocity
% omega(1,1) = 1;
% omega(2,1) = 1;
% the distance between the centers of every two floes when they just touch each other
% the diagonal elements are set to be zero since they represent a floe touch itself
distance_every_two_floes = ones(L,1) * radius + radius' * ones(1,L) - 2 * diag(radius);
% computing the bigger and the smaller radii of each floe pair
radius_max = max(ones(L,1) * radius, radius' * ones(1,L));
radius_min = ones(L,1) * radius + radius' * ones(1,L) - radius_max;
h = thickness';
thickness_min = min(ones(L,1) * thickness, thickness' * ones(1,L));
% Elj and Glj are the Young's modulus (normal direction) and shear modulus
% (tangential direction);
Elj = 50 * ones(L,L);
Glj = 50 * ones(L,L);
% mulj is the contraint coefficient for the tangential force |f_t^lj| <= mulj |f_n^lj| 
mulj = 0.2 * ones(L,L);
% m is the mass 
m = pi * (radius').^2 .* h;     m_truth = m;
% save the contact force
save_contact_force_x = zeros(L,N);
save_contact_force_y = zeros(L,N);
save_contact_force_x_normal = zeros(L,N);
save_contact_force_y_normal = zeros(L,N);
save_contact_force_x_tangential = zeros(L,N);
save_contact_force_y_tangential = zeros(L,N);
save_rotation_force = zeros(L,N);

% ocean drag coefficient, density, portion inside ocean
d_o = 1; rho_o = 1; c_o = 0.9;
alpha_l = d_o * rho_o * pi * radius'.^2;
alpha_L = diag([alpha_l;alpha_l]);
% moment of inertia
I = m .* (radius').^2;
% ocean induced vorticity coefficient
beta_l = d_o * rho_o * pi * (radius').^2 .* (radius').^2;

u_save = zeros(L,N);
v_save = zeros(L,N);
for i = 2:N
    if mod(i,1000) == 0
        disp(i*dt)
    end
    x_loc = [x(:,i-1),y(:,i-1)];

    x(:,i) = x(:,i-1) + (vc_x(:,i-1) + vo_x(:,i-1)) * dt + sqrt(dt) * sigma_x * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + (vc_y(:,i-1) + vo_y(:,i-1)) * dt + sqrt(dt) * sigma_x * randn(L,1); % floe equation in y
    
    vo_x(:,i) = vo_x(:,i-1) + alpha_l ./ m .* (exp(1i * x_loc * kk) * (u_hat(:,i-1) .* transpose(rk(1,:))) - vc_x(:,i-1) - vo_x(:,i-1)) * dt; % velocity induced by the ocean velocity in u
    vo_y(:,i) = vo_y(:,i-1) + alpha_l ./ m .* (exp(1i * x_loc * kk) * (u_hat(:,i-1) .* transpose(rk(2,:))) - vc_y(:,i-1) - vo_y(:,i-1)) * dt; % velocity induced by the ocean velocity in v
    
    u_save(:,i) = exp(1i * x_loc * kk) * (u_hat(:,i-1) .* transpose(rk(1,:)));
    v_save(:,i) = exp(1i * x_loc * kk) * (u_hat(:,i-1) .* transpose(rk(2,:)));
    
    % computing the distance between the centers of different floes
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
    
    
    % finding the floe pairs that have distance smaller than the summation
    % of the two radii (floes overlap with each other). These floes will have resistence forcing to each other. 
    distance_selected = (distance < distance_every_two_floes);
    projection_x = distance_x ./ (distance + eps * eye(L)) .* distance_selected; 
    projection_y = distance_y ./ (distance + eps * eye(L)) .* distance_selected;
    
    % computing the transverse_area (namely the chord length)
    transverse_area = 1 ./ (distance + eps * eye(L)) .* sqrt(4 * distance.^2 .* radius_max.^2 - (distance.^2 - radius_min.^2 + radius_max.^2).^2) .* distance_selected;
    
    
    
    % computing the resistence forcing (normal direction)
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
    normal_force_x = - Elj .* thickness_min .* transverse_area .* (distance_vector_x1_temp + distance_vector_x2_temp + distance_vector_x3_temp); % resistence forcing in the x direction
%     pause

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
    
    
    % computing the normal unit vector (outside direction, pointing towards the other floe)
    normal_direction_x = ((x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)') .* equal_temp_x1 +...
        (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' + 2*pi) .* equal_temp_x2 +...
        (x(:,i-1) * ones(1,L) - ones(L,1) * x(:,i-1)' - 2*pi) .* equal_temp_x3) ./ (distance + eps * eye(L));
    normal_direction_y = ((y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)') .* equal_temp_y1 +...
        (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' + 2*pi) .* equal_temp_y2 +...
        (y(:,i-1) * ones(1,L) - ones(L,1) * y(:,i-1)' - 2*pi) .* equal_temp_y3) ./ (distance + eps * eye(L));
    
    % computing the tangential unit vector (a 90degree rotation from the normal vector; 
    % the positive direction of the tangential vector doesn't matter since the projections
    % to x and y axes will be applied)
    tangential_direction_x = - normal_direction_y;
    tangential_direction_y =   normal_direction_x;
        
    % computing the resistence forcing in the tangential direction
    % computing the velocity difference at the contact interface
    v_diff_x = + (vc_x(:,i-1) + vo_x(:,i-1)) * ones(1,L) - ones(L,1) * (vc_x(:,i-1)' + vo_x(:,i-1)');
    v_diff_y = + (vc_y(:,i-1) + vo_y(:,i-1)) * ones(1,L) - ones(L,1) * (vc_y(:,i-1)' + vo_y(:,i-1)');

    % computing difference between the angular velocity multiplying the radius at the contact interface
    angular_velocity_diff_x = - (ones(L,1) * radius)' .* normal_direction_y' .* (ones(L,1) * omega(:,i-1)')' + (ones(L,1) * radius) .* normal_direction_y .* (ones(L,1) * omega(:,i-1)');
    angular_velocity_diff_y =   (ones(L,1) * radius)' .* normal_direction_x' .* (ones(L,1) * omega(:,i-1)')' - (ones(L,1) * radius) .* normal_direction_x .* (ones(L,1) * omega(:,i-1)');
    
 
    % computing the total tangential force and its projections to the x and
    % y axes
    tangential_force = ((v_diff_x + angular_velocity_diff_x) .* tangential_direction_x + (v_diff_y + angular_velocity_diff_y) .* tangential_direction_y);
    tangential_force_temp = mulj .* sqrt(normal_force_x.^2 + normal_force_y.^2);
    tangential_force_max = min(tangential_force_temp, abs(tangential_force));
    
    tangential_force = ((tangential_force_max == tangential_force_temp) .* tangential_force_temp .* ones(L) + (tangential_force_max == abs(tangential_force)) .* abs(tangential_force) .* ones(L) ).*sign(tangential_force);
    tangential_force_x = Glj .* thickness_min .* transverse_area .* tangential_force .* tangential_direction_x;
    tangential_force_y = Glj .* thickness_min .* transverse_area .* tangential_force .* tangential_direction_y;
    
    % velocity induced by the contact forces
    %save_contact_force_x(:,i-1) = 1./m .* sum(normal_force_x)' + 1./m .* sum(tangential_force_x)';
    %save_contact_force_y(:,i-1) = 1./m .* sum(normal_force_y)' + 1./m .* sum(tangential_force_y)';
    %save_contact_force_x_normal(:,i-1) = 1./m .* sum(normal_force_x)';
    %save_contact_force_y_normal(:,i-1) = 1./m .* sum(normal_force_y)';
    %save_contact_force_x_tangential(:,i-1) = 1./m .* sum(tangential_force_x)';
    %save_contact_force_y_tangential(:,i-1) = 1./m .* sum(tangential_force_y)';

    vc_x(:,i) = 0;% vc_x(:,i-1) + save_contact_force_x(:,i-1) * dt;
    vc_y(:,i) = 0;% vc_y(:,i-1) + save_contact_force_y(:,i-1) * dt;
    % rotation
    t_o = beta_l .* ( exp(1i * x_loc * kk) * ( u_hat(:,i-1) .* transpose( 1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:) ) )/2 - omega(:,i-1) ); 
% t_o
%     pause

    save_rotation_force(:,i-1) = 1./I .* (t_o);% ( sum( (ones(L,1) * radius)' .* (normal_direction_x .* tangential_force_y - normal_direction_y .* tangential_force_x))'  + t_o );
    omega(:,i) = omega(:,i-1) + save_rotation_force(:,i-1) * dt;
    % Periodic boundary conditions
    x(:,i) = mod(x(:,i),2*pi);
    y(:,i) = mod(y(:,i),2*pi);
%     pause
end
v_total_x = vc_x + vo_x;
v_total_y = vc_y + vo_y;
    
% figure 
% hold on
% for l = 1:L
%     th = 0:pi/50:2*pi;
%     xunit = radius(l) * cos(th) + x(l,end);
%     yunit = radius(l) * sin(th) + y(l,end);
%     h = plot(xunit, yunit,'color',[3*radius(l),0.5,0.5]);
% end
% xlim([0, 2*pi ])
% ylim([0, 2*pi ])
% box on

[xx,yy] = meshgrid(linspace(0,2*pi,Dim_Grid), linspace(0,2*pi,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
figure(3) 
for i = 1:100
%     clf
    plot(0,0)
    hold on
    for l = 1:L
        th = 0:pi/50:2*pi;
        xunit = radius(l) * cos(th) + x(l,1+100*(i-1));
        yunit = radius(l) * sin(th) + y(l,1+100*(i-1));
        hh = plot(xunit, yunit,'color',[.2*radius(l),0.5,0.5]);
        text(x(l,1+100*(i-1)),y(l,1+100*(i-1)),num2str(omega(l,1+100*(i-1))));
%         text(x(l,1+100*(i-1)),y(l,1+100*(i-1)),num2str(thickness(l)));
    end
    xlim([0, 2*pi ])
    ylim([0, 2*pi ])
    box on    
    title(['t = ', num2str(dt*100*(i-1))])
    u = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    pause(0.1);
    hold off
    
end