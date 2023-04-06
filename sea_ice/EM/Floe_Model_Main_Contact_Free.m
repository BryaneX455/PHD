% Floe model main file
rng(11);
T = 10; % total time length
dt = 0.001; % time step
N = T/dt; % total time points in time
x = zeros(L,N);
y = zeros(L,N);
x(:,1) = Location(1,:)'; % Lagrangian tracers (floe center location x) evolving with time
y(:,1) = Location(2,:)'; % Lagrangian tracers (floe center location y) evolving with time
vo_x = zeros(L,N); % velocity component x corresponding to the ocean forces
vo_y = zeros(L,N); % velocity component y corresponding to the ocean forces
sigma_x = 0.1/2; % random noise in the tracer equation
omega = zeros(L,N); % angular velocity

% m is the mass 
h = thickness';
m = 1000e9 * pi * (radius').^2 .* h;     m_truth = m;
save_rotation_force = zeros(L,N);

% ocean drag coefficient, density, portion inside ocean
d_o = 5.5; rho_o = 1000e9; c_o = 0.9;
alpha_l = d_o * rho_o * pi * radius'.^2;
alpha_L = diag([alpha_l;alpha_l]);
% moment of inertia
I = m .* (radius').^2;
% ocean induced vorticity coefficient
beta_l = d_o * rho_o * pi * (radius').^2 .* (radius').^2;

for i = 2:N
    if mod(i,1000) == 0
        disp(i*dt)
    end
    x_loc = [x(:,i-1),y(:,i-1)];
    
    % Euler Maruyama
    x(:,i) = x(:,i-1) + (vo_x(:,i-1)) * dt + sqrt(dt) * sigma_x * randn(L,1); % floe equation in x
    y(:,i) = y(:,i-1) + (vo_y(:,i-1)) * dt + sqrt(dt) * sigma_x * randn(L,1); % floe equation in y
    
    vo_x(:,i) = vo_x(:,i-1) + alpha_l ./ m .* 8.64.* (exp(1i * x_loc * kk .*50 /(2*pi)) * (u_hat(:,i-1) .* transpose(rk(1,:))) - vo_x(:,i-1)) * dt; % velocity induced by the ocean velocity in u
    %vo_x(:,i) = vo_x(:,i-1) + alpha_l ./ m .* 8.64.* (exp(1i * x_loc * kk /50 *(2*pi)) * (u_hat(:,i-1) .* transpose(rk(1,:))) - vo_x(:,i-1)) * dt;
    vo_y(:,i) = vo_y(:,i-1) + alpha_l ./ m .* 8.64.* (exp(1i * x_loc * kk .*50 /(2*pi)) * (u_hat(:,i-1) .* transpose(rk(2,:))) - vo_y(:,i-1)) * dt; % velocity induced by the ocean velocity in v
    %vo_y(:,i) = vo_y(:,i-1) + alpha_l ./ m .* 8.64.* (exp(1i * x_loc * kk /50 *(2*pi)) * (u_hat(:,i-1) .* transpose(rk(2,:))) - vo_y(:,i-1)) * dt;
   
    % rotation
    t_o = beta_l .* ( exp(1i * x_loc * kk *50/ 2/ pi) * ( u_hat(:,i-1) .* transpose( 1i * rk(2,:) .* kk(2,:) - 1i * rk(1,:) .* kk(1,:) ) )/2 - omega(:,i-1) ); 
    save_rotation_force(:,i-1) = 1./I .* (t_o); % ( sum( (ones(L,1) * radius)' .* (normal_direction_x .* tangential_force_y - normal_direction_y .* tangential_force_x))'  + t_o );
    omega(:,i) = omega(:,i-1) + save_rotation_force(:,i-1) * dt;
    % Periodic boundary conditions
    x(:,i) = mod(x(:,i),50);
    y(:,i) = mod(y(:,i),50);
% pause
end
v_total_x = vo_x;
v_total_y = vo_y;

[xx,yy] = meshgrid(linspace(0,50,Dim_Grid), linspace(0,50,Dim_Grid));
x_vec = [reshape(xx,[],1), reshape(yy,[],1)]; 
figure 
for i = 1:100
    plot(0,0)
    hold on
    for l = 1:L
        th = 0:pi/50:2*pi;
        xunit = radius(l) * cos(th) + x(l,1+100*(i-1));
        yunit = radius(l) * sin(th) + y(l,1+100*(i-1));
        hh = plot(xunit, yunit,'color',[.2*radius(l),0.5,0.5]);
        text(x(l,1+100*(i-1)),y(l,1+100*(i-1)),num2str(omega(l,1+100*(i-1))));
    end
    xlim([0, 50])
    ylim([0, 50])
    box on    
    title(['t = ', num2str(dt*100*(i-1))])
    u = exp(1i * x_vec * kk /2/pi*50) * (u_hat(:,1+100*(i-1)) .* transpose(rk(1,:)));
    v = exp(1i * x_vec * kk /2/pi*50) * (u_hat(:,1+100*(i-1)) .* transpose(rk(2,:)));
    u = reshape(u, Dim_Grid, Dim_Grid);
    v = reshape(v, Dim_Grid, Dim_Grid);
    quiver(xx, yy, u, v, 'linewidth',1)
    pause(0.1);
    hold off    
end