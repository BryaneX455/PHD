% Chapter 9
% Parameter estimation with partial observations via expectation-maximization
% algorithm; the test model here is the noisy Lorenz 84 model, where the
% observed variables are y and z while x is unobserved

rng(11) % fix the random number seed to reproduce the results
N = 10000; % total number of numerical integration time step
dt = 0.001; % numerical integration time step
% state variables
x_truth = zeros(1,N);
y_truth = zeros(1,N);
z_truth = zeros(1,N);
% model parameters
g_truth = 1;
b_truth = 4;
a_truth = 1/4;
f_truth = 2;
% noise coefficients (two regimes with values = .1 or 3)
sigma_x_truth = .1;
sigma_y_truth = .1;
sigma_z_truth = .1;

% sigma_x_truth = 3;
% sigma_y_truth = 3;
% sigma_z_truth = 3;

% generating the true signal
for i = 2:N
    x_truth(i) = x_truth(i-1) + ( - (y_truth(i-1)^2 + z_truth(i-1)^2) - a_truth * x_truth(i-1) + f_truth) * dt + sigma_x_truth * randn * sqrt(dt);
    y_truth(i) = y_truth(i-1) + ( - b_truth * x_truth(i-1) * z_truth(i-1) + x_truth(i-1) * y_truth(i-1) - y_truth(i-1) + g_truth) * dt + sigma_y_truth * randn * sqrt(dt);
    z_truth(i) = z_truth(i-1) + ( b_truth * x_truth(i-1) * y_truth(i-1) + x_truth(i-1) * z_truth(i-1) - z_truth(i-1)) * dt + sigma_z_truth * randn * sqrt(dt);
end
% showing the true signal and the associated statistics
% The short time series of y and z here are utilized for parameter 
% estimation, while at the end of this code long time series are utilized
% for the comparison of the model behiavor between the models with the 
% true parameters and the estimated parameters
figure
subplot(3,1,1)
plot(dt:dt:N*dt,x_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(a) The associated unobserved x','fontsize',12)
subplot(3,1,2)
plot(dt:dt:N*dt,y_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(b) Observed y for parameter estimation','fontsize',12)
subplot(3,1,3)
plot(dt:dt:N*dt,z_truth,'b','linewidth',2)
box on
set(gca,'fontsize',12)
title('(c) Observed z for parameter estimation','fontsize',12)

%%%%%%%%%%%%%%%%%%% Parameter Estimation via EM %%%%%%%%%%%%%%%%%%%%%%%%

% initial guess of the parameters
g = 1*2;
b = 4*2;
a = 1/4*2;
f = 2*2;

Theta = [g;b;a;f]; % parameters in the deterministic part
R = diag([1^2,1^2,1^2 ])*dt; % covariance in the EM algorithm, related to the parameters in the stochastic part
sigma_x = sqrt(R(3,3)/dt); % put the noise coefficient of the unobserved variable x in R(3,3)
sigma_y = sqrt(R(1,1)/dt);
sigma_z = sqrt(R(2,2)/dt);
n1 = length(Theta); % number of parameters in the deterministic part
n2 = 3; % number of parameters in the stochastic part
KK = 200; % Total number of EM iterations

Param_save = zeros(n1+n2,KK); % save the parameters

Param_save(:,1) = [Theta;sigma_x;sigma_y;sigma_z]; % including initial guess of the parameters
Param_truth = [g_truth;b_truth;a_truth;f_truth;sigma_x_truth;sigma_y_truth;sigma_z_truth]; % the true parameter values

% EM algorithm
for k = 2:KK
    if mod(k,10) == 1
        disp(['True parameter values'])
        disp('g,   b,   a,   f,   \sigma_x,   \sigma_y,   \sigma_z')
        disp(Param_truth')
        disp(['Parameter values at the ', num2str(k), '-th iteration'])
        disp(Param_save(:,k-1)')
    end
    % E Step; nonlinear smoothing    
    gamma_mean_trace = zeros(1,N); % filtering mean
    gamma_cov_trace = zeros(1,N); % filtering variance
    gamma_cov_trace(:,1) = 1;
    mu_s = zeros(1,N); % smoothing mean
    R_s = zeros(1,N); % smoothing variance
    C = zeros(1,N); % auxiliary matrix in smoothing

    % filtering
    gamma_mean0 = gamma_mean_trace(:,1);
    gamma_cov0 = gamma_cov_trace(:,1);
    for i = 2:N
        % observations: y and z
        u0 = [y_truth(i-1);z_truth(i-1)];
        u = [y_truth(i);z_truth(i)];
        
        % matrices and vectors needed in filtering
        a1 = - a;
        a0 = - y_truth(i-1)^2 - z_truth(i-1)^2 + f;        
        A0 = [- y_truth(i-1) + g;
            - z_truth(i-1)];
        A1 = [- b * z_truth(i-1) + y_truth(i-1);
            b * y_truth(i-1) + z_truth(i-1)];
        b1 = sigma_x;
        invBoB = [1/sigma_y^2, 0;
            0, 1/sigma_z^2];

        % filtering updating
        gamma_mean = gamma_mean0 + (a0 + a1 * gamma_mean0) * dt + (gamma_cov0 * A1') * invBoB * (u-u0 - A0*dt-A1 * gamma_mean0 * dt);
        gamma_cov = gamma_cov0 + (a1 * gamma_cov0 + gamma_cov0 * a1' + b1 * b1' - (gamma_cov0 * A1') * invBoB * (gamma_cov0 * A1')') * dt;     
        % saving the filtering results
        gamma_mean_trace(i) = gamma_mean;
        gamma_cov_trace(i) = gamma_cov;

        gamma_mean0 = gamma_mean;
        gamma_cov0 = gamma_cov;
    end
    % backward smoothing
    mu_s(end) = gamma_mean; % the starting point is the final value of the filtering
    R_s(end) = gamma_cov;
    for i = N-1:-1:1
        % observations
        u0 = [y_truth(i);z_truth(i)];
        % matrices and vectors needed in smoothing
        a0 = - y_truth(i)^2 - z_truth(i)^2 + f;
        b1 = sigma_x;
        C(i) = gamma_cov_trace(i) * (1 + a1 * dt)' * (b1 * b1' * dt + (1 + a1 * dt) * gamma_cov_trace(i) * (1 + a1 * dt)')^(-1);
        mu_s(i) = gamma_mean_trace(i) + C(i) * (mu_s(i+1) - a0 * dt - ( 1 + a1 * dt) * gamma_mean_trace(i)); % smoothing mean
        R_s(i) = gamma_cov_trace(i) + C(i) * (R_s(i+1) - (1 + a1 * dt) * gamma_cov_trace(i) * (1 + a1 * dt)' - b1 * b1 * dt) * C(i)'; % smoothing variance
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

    % M step via a quadratic optimization
    % updating the drift terms
    for i = 1+10:N-1-10
        X = mu_s(i);
        XX = mu_s(i)^2 + R_s(i);
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
    
    % updating diffusion terms
    flag_num = 0;
    R = zeros(n2,n2);
    Term1_total = 0;
    Term2_total = 0;
    Term3_total = 0;
    Term4_total = 0;
    for i = 1+10:N-1-10
        X = mu_s(i);
        XX = mu_s(i)^2 + R_s(i);
        XXj = mu_s(i) * mu_s(i+1) + R_s(i+1) * C(i)';
        a1 = - a;
        a0 = - y_truth(i)^2 - z_truth(i)^2 + f;
        
        A0 = [- y_truth(i) + g;
            - z_truth(i)];
        A1 = [- b * z_truth(i) + y_truth(i);
            b * y_truth(i) + z_truth(i)];
        YZj1 = [y_truth(i+1); z_truth(i+1)];
        YZj = [y_truth(i); z_truth(i)];
        Term1 = [(YZj1 - YZj) * (YZj1 - YZj)', (YZj1 - YZj) * (mu_s(i+1) - mu_s(i))';
            (mu_s(i+1) - mu_s(i)) * (YZj1 - YZj)', mu_s(i+1)^2 + R_s(i+1) + XX - XXj - XXj'];
        T121 = (mu_s(i+1) - mu_s(i)) * A0' + XXj * A1' - XX * A1';
        T122 = (mu_s(i+1) - mu_s(i)) * a0' + XXj * a1' - XX * a1';
        Term2 = [(YZj1 - YZj) * (A0 + A1 * mu_s(i))', (YZj1 - YZj) * (a0 + a1 * mu_s(i))';
            T121, T122] * dt;
        Term3 = Term2';
        T211 = A0 * A0' + A0 * X' * A1' + A1 * X * A0' + A1 * XX * A1';
        T212 = A0 * a0' + A0 * X' * a1' + A1 * X * a0' + A1 * XX * a1';
        T221 = a0 * A0' + a0 * X' * A1' + a1 * X * A0' + a1 * XX * A1';
        T222 = a0 * a0' + a0 * X' * a1' + a1 * X * a0' + a1 * XX * a1';
        Term4 = [T211, T212;
            T221, T222] * dt^2;
        flag_num = flag_num + 1;
        R = R + (Term1 - Term2 - Term3 + Term4);R = diag(diag(R));
        Term1_total = Term1_total + Term1;
        Term2_total = Term2_total + Term2;
        Term3_total = Term3_total + Term3;
        Term4_total = Term4_total + Term4;
    end
    g = Theta(1);
    b = Theta(2);
    a = Theta(3);
    f = Theta(4);
    R = R/flag_num;
    Param_save(:,k) = [Theta;sqrt(R(3,3)/dt);sqrt(R(1,1)/dt);sqrt(R(2,2)/dt)]; 
    sigma_y = Param_save(6,k);
    sigma_z = Param_save(7,k);
    sigma_x = sigma_x_truth; %%% here we fix the noise coefficient in the unobserved process.
    % Unlike the noise coefficients in the observed processes which can be 
    % updated using quadratic variation, the noise in the unobserved 
    % process is hard to be updated directly. Usually this parameter is 
    % updated by using a change of measure. 
    % In the presence of a small noise coefficient, the direct method can 
    % sometimes be used but otherwise the observability issue makes the 
    % convergence very slow.    
%     sigma_x = Param_save(5,k);
    % One method to accelerate the convergence is as follows:
%     eta = 50;
%     sigma_x = Param_save(5,k-1) + eta * (Param_save(5,k) - Param_save(5,k-1));Param_save(5,k)  = sigma_x;

end
% showing the results
figure
for i = 1:n1+n2
    subplot(2,4,i)
    hold on
    plot(1:KK, Param_save(i,:),'b','linewidth',2)
    plot([1,KK], [Param_truth(i), Param_truth(i)],'--k','linewidth',2)
    box on
    set(gca,'fontsize',12)
    if i == 1
        title('g','fontsize',14)
    elseif i == 2
        title('b','fontsize',14)
    elseif i == 3
        title('a','fontsize',14)
    elseif i == 4
        title('f','fontsize',14)
    elseif i == 5
        title('\sigma_x','fontsize',14)
    elseif i == 6
        title('\sigma_y','fontsize',14)
    elseif i == 7
        title('\sigma_z','fontsize',14)
    end
    xlabel('k');
end

% the estimated parameters, which are the final value of the iteration
g_est = Param_save(1,end);
b_est = Param_save(2,end);
a_est = Param_save(3,end);
f_est = Param_save(4,end);
sigma_x_est = Param_save(5,end);
sigma_y_est = Param_save(6,end);
sigma_z_est = Param_save(7,end);
figure
subplot(3,1,1)
hold on
plot(dt:dt:N*dt,x_truth,'b','linewidth',2)
plot(dt:dt:N*dt,mu_s_2,'c','linewidth',2)
patch([dt:dt:N*dt,N*dt:-dt:dt],[mu_s_2+2*sqrt(R_s_2),mu_s_2(end:-1:1)-2*sqrt(R_s_2(end:-1:1))],'c','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',12)
legend('Truth of x','Smoother')
title('After k = 2 iteration steps')
subplot(3,1,2)
hold on
plot(dt:dt:N*dt,x_truth,'b','linewidth',2)
plot(dt:dt:N*dt,mu_s_5,'c','linewidth',2)
patch([dt:dt:N*dt,N*dt:-dt:dt],[mu_s_5+2*sqrt(R_s_5),mu_s_5(end:-1:1)-2*sqrt(R_s_5(end:-1:1))],'c','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',12)
title('After k = 5 iteration steps')
subplot(3,1,3)
hold on
plot(dt:dt:N*dt,x_truth,'b','linewidth',2)
plot(dt:dt:N*dt,mu_s,'c','linewidth',2)
patch([dt:dt:N*dt,N*dt:-dt:dt],[mu_s+2*sqrt(R_s),mu_s(end:-1:1)-2*sqrt(R_s(end:-1:1))],'c','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',12)
title('At the last iteration steps')

N2 = N*50;
x_truth2 = zeros(1,N2);
y_truth2 = zeros(1,N2);
z_truth2 = zeros(1,N2);
x_est = zeros(1,N);
y_est = zeros(1,N);
z_est = zeros(1,N);
% comparing the model trajectories and statistics between the one with the
% true parameters and the one with the estimated parameters
for i = 2:N2
    x_truth2(i) = x_truth2(i-1) + ( - (y_truth2(i-1)^2 + z_truth2(i-1)^2) - a_truth * x_truth2(i-1) + f_truth) * dt + sigma_x_truth * randn * sqrt(dt);
    y_truth2(i) = y_truth2(i-1) + ( - b_truth * x_truth2(i-1) * z_truth2(i-1) + x_truth2(i-1) * y_truth2(i-1) - y_truth2(i-1) + g_truth) * dt + sigma_y_truth * randn * sqrt(dt);
    z_truth2(i) = z_truth2(i-1) + ( b_truth * x_truth2(i-1) * y_truth2(i-1) + x_truth2(i-1) * z_truth2(i-1) - z_truth2(i-1)) * dt + sigma_z_truth * randn * sqrt(dt);
end

for i = 2:N2
    x_est(i) = x_est(i-1) + ( - (y_est(i-1)^2 + z_est(i-1)^2) - a_est * x_est(i-1) + f_est) * dt + sigma_x_est * randn * sqrt(dt);
    y_est(i) = y_est(i-1) + ( - b_est * x_est(i-1) * z_est(i-1) + x_est(i-1) * y_est(i-1) - y_est(i-1) + g_est) * dt + sigma_y_est * randn * sqrt(dt);
    z_est(i) = z_est(i-1) + ( b_est * x_est(i-1) * y_est(i-1) + x_est(i-1) * z_est(i-1) - z_est(i-1)) * dt + sigma_z_est * randn * sqrt(dt);
end
Lag = 5000; % lag in computing the ACF
figure
for i = 1:3
    if i == 1
        variable1 = x_truth2;
        variable2 = x_est;
    elseif i == 2
        variable1 = y_truth2;
        variable2 = y_est;
    elseif i == 3
        variable1 = z_truth2;
        variable2 = z_est;
    end
    subplot(3,5,[1:3]+5*(i-1))
    hold on
    plot(dt:dt:N2*dt,variable1,'b','linewidth',2)
    plot(dt:dt:N2*dt,variable2,'r','linewidth',2)
    box on
    set(gca,'fontsize',16)
    if i == 1
        title('Time series of x')
    elseif i == 2
        title('Time series of y')
    else
        title('Time series of z')
        xlabel('t')
    end
    subplot(3,5,4+5*(i-1))
    ACF1 = autocorr(variable1, Lag);
    ACF2 = autocorr(variable2, Lag);
    hold on
    plot(0:dt:Lag*dt, ACF1,'b','linewidth',2)
    plot(0:dt:Lag*dt, ACF2,'r','linewidth',2)
    box on
    set(gca,'fontsize',16)
    if i == 1
        title('ACF')
    end
    if i == 3
        xlabel('t')
    end
    subplot(3,5,5+5*(i-1))
    [fi1,xx] = ksdensity(variable1(1:10:end));
    [fi2,xx] = ksdensity(variable2(1:10:end),xx);
    hold on
    plot(xx,fi1,'b','linewidth',2)
    plot(xx,fi2,'r','linewidth',2)
    box on
    set(gca,'fontsize',16)
    if i == 1
        title('PDF')
    end
end