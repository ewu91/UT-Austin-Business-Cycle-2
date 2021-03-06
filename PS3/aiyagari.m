%{
Problem Set 3 - Business Cycle Analysis - Spring 2021 - UT Austin
Prof: Chris Boehm
Edson An An Wu
Date: Mar, 1st, 2021

Replication of Aiyagari, with transition matrix being:
[1/2 1/2 0; 1/4 1/2 1/4; 0 1/2 1/2]
%}
clear all
close all
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

beta = 0.95;
alpha = 0.36;
delta = 0.1;

%--------------------------------------------------------------------------
% Income Process
%--------------------------------------------------------------------------
ygrid = [0.5, 1.0, 1.5];
ny = length(ygrid);
P =[1/2 1/2 0; 1/4 1/2 1/4; 0 1/2 1/2];

%--------------------------------------------------------------------------
% Asset Grid and etc
%--------------------------------------------------------------------------
amin = 0;
amax = 100;
agrid = amin:0.05:amax;
%agrid = amin:1:amax;
na = length(agrid);

%--------------------------------------------------------------------------
% Pre-Allocating Variables
%--------------------------------------------------------------------------
g = zeros(na,ny); % savings policy function
gc = zeros(na,ny); % consumption policy function
v0 = zeros(na,ny);
v1 = zeros(na,ny); % value function
u = zeros(na,ny,na);
rmin = 0.01;
rmax = 0.05;
L = 1;

T = 11000;
g_timeseries = zeros(T,1);
y_timeseries = zeros(T,1);
g_timeseries(1) = 3;
y_timeseries(1) = 1.0;

%--------------------------------------------------------------------------
% Technical Stuff
%--------------------------------------------------------------------------
tol = 0.0001;
maxit = 10^5;
iter = 0;
error = tol + 1;
iter_r = 0;
error_r = tol + 1;
%%
tic
% Market Clearing: Loop to calculate r equilibrium
while (abs(error_r) > tol) && (iter_r < maxit)
    
    % Updating r0
    r0 = 0.5*rmin + 0.5*rmax;    
    disp(['iteration for r: ', num2str(iter_r), ', error: ', num2str(error_r), ', r value:', num2str(r0), ', rmax: ', num2str(rmax), ', rmin:', num2str(rmin)])        
    R0 = r0 + delta; % Rental rate
    K0 = ((R0/alpha)*(1/L^(1-alpha)))^(1/(alpha-1)); % Capital    
    w0 = (1-alpha)*K0^alpha*L^(-alpha); % Wage

    % Calculating flow utility for all possible combinations of (a,y,a')    
    for ia = 1:na
        for iy = 1:ny            
            c = w0*ygrid(iy) + agrid(ia) - agrid(:)/(1+r0);
            u(ia,iy,c<=0) = - Inf;
            u(ia,iy,c>0) = log(c(c>0));                                            
        end
    end            
    
    error = tol+1;
    iter = 0;
    % Fixed Point for the value function
    while (error > tol) && (iter < maxit)
        % Loop over income
        for iy = 1:ny            

            % Calculating the expected value for each value of a'
            Ev0 = P(iy,:)*v0';

            % Calculating the current utility for each value of a'
            u_temp = reshape(u(:, iy, :),[na,na]);

            % Calculating the value function
            vtemp = u_temp + repmat(beta*Ev0,na,1);

            % Calculating the v1
            [v1(:,iy), maxindex] = max(vtemp,[],2);
            g(:,iy) = agrid(maxindex);

        end
        error = max(abs(v1(:)-v0(:)));
        iter = iter + 1;
        if mod(iter,25) == 0
            disp(['iteration: ', num2str(iter), ', error: ', num2str(error)])        
        end
        v0 = v1;

    end
    
    for iy = 1:ny
        gc(:,iy) = w0*ygrid(iy) + agrid(:) - g(:,iy)/(1+r0);
    end
    
    rng(123);
    % Simulating Time Series of policy functions to get aggregate capital
    for t=2:T
        current_a = g_timeseries(t-1);
        current_y = y_timeseries(t-1);
        iy = find(ygrid == current_y);
        next_a = interp1(agrid,g(:,iy)',current_a);
        
        draw_y = rand(1,1);        
        iy = 1*(draw_y < P(iy,1)) + 2*(draw_y < P(iy,1) + P(iy,2))*(draw_y > P(iy,1)) + 3*(draw_y > P(iy,1) + P(iy,2));
        g_timeseries(t) = next_a;
        y_timeseries(t) = ygrid(iy);    
    end
    % Updating Aggregate Capital
    K1 = mean(g_timeseries(1001:T));
    R1 = alpha*K1^(alpha-1)*L^(1-alpha);
    r1 = R1 - delta; % new interest rate
    w1 = (1-alpha)*K1^alpha*L^(-alpha);
    
    iter_r = iter_r + 1;
    % Bissection part
    error_r = r1 - r0;
    if error_r > 0
        rmin = r0;
    else
        rmax = r0;
    end
    
end    
toc    
%%
%==========================================================================
% Plotting Figures
%==========================================================================

figure(1)
set(gcf,'Position',[100 100 1200 720])
hold on
plot(agrid,g(:,1),'DisplayName','y=0.5')
plot(agrid,g(:,2),'DisplayName','y=1')
plot(agrid,g(:,3),'DisplayName','y=1.5')
hold off
title('Policy Function for Assets')
xlabel('Asset State')
ylabel('Next Period Assets')
legend('Location','southeast')
orient(gcf,'landscape')
saveas(gcf, 'aiyagari_policy_function.pdf')

figure(2)
set(gcf,'Position',[100 100 1200 720])
hold on
plot(agrid,gc(:,1),'DisplayName','y=0.5')
plot(agrid,gc(:,2),'DisplayName','y=1')
plot(agrid,gc(:,3),'DisplayName','y=1.5')
hold off
title('Policy Function for Consumption')
xlabel('Asset State')
ylabel('Consumption')
legend('Location','southeast')
orient(gcf,'landscape')
saveas(gcf, 'aiyagari_consumption.pdf')

figure(3)
set(gcf,'Position',[100 100 1200 720])
hold on
plot(agrid,v1(:,1),'DisplayName','y=0.5')
plot(agrid,v1(:,2),'DisplayName','y=1')
plot(agrid,v1(:,3),'DisplayName','y=1.5')
hold off
title('Value Function')
xlabel('Asset State')
ylabel('Value Function')
legend('Location','southeast')
orient(gcf,'landscape')
saveas(gcf, 'aiyagari_value_function.pdf')

figure(5)
set(gcf,'Position',[100 100 1200 720])
hold on
histogram(g_timeseries(1001:T))
hold off
title('Stationary Distribution of Assets')
xlabel('Asset State')
ylabel('Frequency')
orient(gcf,'landscape')
saveas(gcf, 'aiyagari_value_function.pdf')
    
    
    
    









