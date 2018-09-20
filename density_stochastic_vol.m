% Script for computing the p.d.f. of stochastic volatility model

clear;

% parmaeter setting
kappa = 11; theta = 0.2; gam = 1.5;
a = 1;      % possible values for a = 0, 1
b = 1;      % possible values for b = 1/2, 1, 3/2


% for numerical computation
dt = 1/5/250;  dx = 0.0025; dy = 0.001;
tol_x = 1.0000e-5;  tol_y = 1.0000e-8; % tolerance for transition probability density function

N = 50; T = N*dt;

% setting for x domain
x = (0.2*theta+dx):dx:1.8*theta;    % possible presetting for domain for x, heuristic

% dynamic setting for y domain and corresponding f
Num_y = 250;

y = 0;
f = normpdf(y + x, x+ kappa*x.^a.*( theta - x)*dt, gam*x.^b*sqrt(dt))';
f_test_right = f; f_test_left = f;

while (max(f_test_right) > tol_y || max(f_test_left) > tol_y) || length(y) < Num_y
    y = [min(y) - dy, y, max(y) + dy];
    f_test_right = normpdf(max(y) + x, x+ kappa*x.^a.*( theta - x )*dt, gam*x.^b*sqrt(dt))';
    f_test_left  = normpdf(min(y) + x, x+ kappa*x.^a.*( theta - x )*dt, gam*x.^b*sqrt(dt))';
    f = [f_test_left, f, f_test_right];
end

% if the number of steps in y is too larger, reduce it.
[y, dy, f] = redefine_yf(x, y, f, Num_y, 'spline');

% set appropriate num_chi
center_chi = theta;                    % this is arbitrarily assumed.
min_chi = center_chi; max_chi = center_chi;
num_chi = 1;
while normpdf(min_chi, center_chi + kappa*center_chi^a*( theta - center_chi)*dt, gam*max(x)^b*sqrt(dt)) > tol_x || ...
        normpdf(max_chi, center_chi + kappa*center_chi^a*( theta - center_chi)*dt, gam*max(x)^b*sqrt(dt)) > tol_x
    min_chi = min_chi - dx; max_chi = max_chi + dx; num_chi = num_chi + 2;
end

% Pre calculation
chi = zeros(num_chi, length(x));
pdfs = zeros(num_chi, length(x));
chi_idx_xs = zeros(num_chi, length(x));
for i = 1:length(x)
    % set a vector chi for possible range of X_{n+1} around X_n = x(i) 
    chi(:, i) = x(i) - floor(num_chi/2) * dx : dx : x(i) + floor(num_chi/2) * dx;
    pdfs(:, i) = normpdf(chi(:,i), x(i) + kappa*x(i).^a.*( theta - x(i))*dt, gam*x(i).^b*sqrt(dt) );
    chi_idx_xs(:, i) = simple_query(chi(:,i), x);
end


% main loop over time
for n = N-1:-1:1
    new_f = zeros(length(x), length(y));
    for i=1:length(x)
        % main part
        for j=1:length(y)
            y_minus_h = y(j) - (chi(:, i) - x(i));
            integrand_psi = f(sub2ind(size(f), chi_idx_xs(:, i), simple_query(y_minus_h, y)));
            new_f(i,j) = sum(integrand_psi.*pdfs(:, i))*dx;   %simple integration works fine
        end
    end

    % Extend grid, if the side of pdf is not enough close to zero
    original_y = y;
    while (max(new_f(:,end)) > tol_y || max(new_f(:,1)) > tol_y)  
        new_j = [];
        if max( abs(new_f(:,1) )) > tol_y
            y = [min(y) - dy, y];
            new_f = [zeros(length(x),1), new_f];
            new_j = 1;
        end
        if max( abs(new_f(:,end) )) > tol_y
            y = [y, max(y) + dy];
            new_f = [new_f, zeros(length(x),1)];
            new_j = [new_j, length(y)];
        end
        for i=1:length(x)
            % main part repeated
            for j_=1:length(new_j)
                j = new_j(j_);
                y_minus_h = y(j) - (chi(:,i) - x(i));
                integrand_psi = f(sub2ind(size(f), chi_idx_xs(:,i), simple_query(y_minus_h, original_y))); 
                new_f(i,j) = sum(integrand_psi.*pdfs(:, i))*dx;
            end
        end
    end
    
    f = new_f;
    % redefine grid if Num_y is too large
    [y, dy, f] = redefine_yf(x, y, f, Num_y, 'spline');
    
end

% choose pdf
V0 = theta;
pdf = f(abs(x-theta)<dx/2,1:1:end);
% numerical part finished


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

