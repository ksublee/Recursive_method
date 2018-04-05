% Script for computing the p.d.f. of stochastic volatility model

clear;

% parmaeter setting
kappa = 11; theta = 0.2; gam = 1.5;

a = 1;      % possible value for a = 0, 1
b = 1;      % possible value for b = 1/2, 1, 3/2

dt = 1/5/250;  dx = 0.0025; dy = 0.001;

tol_x = 1.0000e-5;  tol_y = 1.0000e-8; % tolerance for transition probability density function

N = 50;
T = N*dt;

% setting for x domain
x = (0.2*theta+dx):dx:1.8*theta;    % possible presetting for domain for x


% xx is possible range of X_{N} with X_{N-1} = x_start
x_start = theta;                    % this is arbitrarily assumed.
xx = theta;                         
while normpdf(min(xx), x_start + kappa*x_start^a*( theta - x_start)*dt, gam*max(x)^b*sqrt(dt)) > tol_x || ...
        normpdf(max(xx), x_start + kappa*x_start^a*( theta - x_start)*dt, gam*max(x)^b*sqrt(dt)) > tol_x
    xx = [min(xx)-dx, xx, max(xx)+dx];
end


% dynamic setting for y domain

Num_y = 250;

% first pdf grid, f
y_value = 0;
f_vector = normpdf(y_value + x, x+ kappa*x.^a.*( theta - x)*dt, gam*x.^b*sqrt(dt));
f_test_right = f_vector; f_test_left = f_vector;
f = f_vector';
y = y_value;

while (max(f_test_right) > tol_y || max(f_test_left) > tol_y) || length(y) < Num_y
    y = [min(y) - dy, y, max(y) + dy];
    f_test_right = normpdf(max(y) + x, x+ kappa*x.^a.*( theta - x )*dt, gam*x.^b*sqrt(dt));
    f_test_left  = normpdf(min(y) + x, x+ kappa*x.^a.*( theta - x )*dt, gam*x.^b*sqrt(dt));
    f = [f_test_left', f, f_test_right'];
end

% if the number of steps in y is too larger, reduce it.
new_y = y;
if length(y) > Num_y+ 1
    dy = (max(y)-min(y))/Num_y;
    new_y = min(y):dy:max(y);
    [mesh_y, mesh_x] = meshgrid(new_y, x);
    f = interp2(y, x, f, mesh_y, mesh_x, 'spline');
    y = new_y;
end

% loop for time
for n = N-1:-1:1

    new_f = zeros(length(x), length(y));
    min_y = min(y);
    
    for i=1:length(x)
        
        % possible range of X_{n+1} with X_n = x(i) 
        x1 = (x(i) - xx(round(length(xx)/2))) + xx;
        
        % if query points are out of bound, take maximum or minimum
        q_points_x = x1;
        q_points_x(q_points_x > max(x)) = max(x);
        q_points_x(q_points_x < min(x)) = min(x);
        q_x_index = round((q_points_x - min(x)) / dx) + 1;  % index

        % approximated pdf over x1
        transition_pdf = normpdf(x1, x(i) + kappa*x(i).^a.*( theta - x(i))*dt, gam*x(i).^b*sqrt(dt) );
        
        for j=1:length(y)
            
            q_points = y(j) - (x1 - x(i));
            q_y_index = round((q_points - min_y)/dy)+1;
            
            integrand_f = zeros(1, length(q_points));
            
            for m=1:length(q_points)
                if q_y_index(m) >= 1 && q_y_index(m) <= Num_y
                    integrand_f(m) = f(q_x_index(m), q_y_index(m));  % reference from previous f
                end
            end                                   
            new_f(i,j) = sum(integrand_f.*transition_pdf)*dx;   %simple integration works fine
        end
    end

    % Extend grid, if the pdf is not enough close to zero
    original_y = y;
    while (max(new_f(:,end)) > tol_y || max(new_f(:,1)) > tol_y)  
        j_vector = [];
        if max( abs(new_f(:,1) )) > tol_y
            y = [min(y) - dy, y];
            new_f = [zeros(length(x),1), new_f];
            j_vector = [1,j_vector];
        end
        if max(abs(new_f(:,end) )) > tol_y
            y = [y, max(y) + dy];
            new_f = [new_f, zeros(length(x),1)];
            j_vector = [j_vector, length(y)];
        end

        for i=1:length(x)
            x1 = (x(i) - xx(round(length(xx)/2))) + xx;

            q_points_x = x1;
            q_points_x(q_points_x > max(x)) = max(x);
            q_points_x(q_points_x < min(x)) = min(x);
            q_x_index = round((q_points_x - min(x))/dx)+1;
            
            transition_pdf = normpdf(x1, x(i) + kappa*x(i).^a.*( theta - x(i))*dt, gam*x(i).^b*sqrt(dt) );
            for j_=1:length(j_vector)
                j=j_vector(j_);
                q_points = y(j) - (x1 - x(i));
                q_y_index = round((q_points - min(original_y))/dy)+1;
            
                integrand_f = zeros(1, length(q_points));
                for m=1:length(q_points)
                    if q_y_index(m) >= 1 && q_y_index(m) <= length(original_y)
                        integrand_f(m) = f(q_x_index(m), q_y_index(m));
                    end
                end               
                new_f(i,j) = trapz(x1, integrand_f.*transition_pdf);
            end
        end
    end
    
    f = new_f;
    
    new_y = y;
    if length(y) > Num_y + 1
        dy = (max(y)-min(y))/Num_y;
        new_y = min(y):dy:max(y);
        [mesh_y, mesh_x] = meshgrid(new_y, x);
        f = interp2(y, x, f, mesh_y, mesh_x, 'spline');
        y = new_y;
    end
    
end

% choose pdf
V0 = theta;
pdf = f(abs(x-theta)<dx/2,1:1:end);



%simulation

Num_simul = 100000;
V = zeros(1, N+1);
V(1) = V0;
last_V = zeros(1, Num_simul);
ds = 1/5/250;
Ns = round(T/ds);

for j=1:Num_simul
    W = normrnd(0,1,1,Ns);    
    for i = 2:Ns+1
        V(i) = V(i-1) + kappa*V(i-1).^a.*(theta-V(i-1))*ds + gam*V(i-1)^b*W(i-1)*sqrt(ds);
    end
    last_V(j) = V(end)-V0;
end


new_y = min(real(last_V))-10*dy:dy:max(real(last_V));
new_pdf = interp1(y, pdf, new_y, 'pchip', 0);
new_cdf = cumsum(new_pdf)*dy;


figure(2)
plot(new_y+0.2, new_pdf, 'LineWidth', 1)
hold on
histogram(real(last_V)+0.2, 40, 'Normalization', 'pdf', 'EdgeColor', 'k', 'FaceColor', 'none', 'LineStyle','-')
hold off

