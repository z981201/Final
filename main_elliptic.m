clc
clear
close all
addpath(genpath('..'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Initialize
% 0.1 Gamma
a = 2; % major axis
b = 1; % minor axis
gamma = {@(t,r) a*r.*cos(t);
         @(t,r) b*r.*sin(t)}; 
gamma_dt = {@(t,r) -r.*a*sin(t);
            @(t,r)  r.*b*cos(t)}; 
gamma_dt2 = {@(t) -1.*a*cos(t);
            @(t)  -1.*b*sin(t)}; 

% 0.2 u_exact 
u_exact = rand(3,1);
u_exact = u_exact/sum(u_exact);
u_exact = @(x,y) 5+u_exact(1)/(a*b)*x.*y +...
                 u_exact(2)/exp(b)*exp(-y).*sin(x) + u_exact(3)/exp(a)*exp(x).*cos(y);

% 0.3 src in Gamma
N = 100; % num of src
theta = (0:N-1)'*2*pi/N;
h = 2*pi/N;
src = [gamma{1}(theta,1) gamma{2}(theta,1)];

% 0.4 target in Omega
m = 400; % num of d_{rho}
rho_correct = 0.5; % correct by power function
rho = ((1:m)/(m+1)).^(1/(1+rho_correct));
% rho_correct = (N-1)/N; % rand and correct by linear
% rho = rand(1,m)*rho_correct;
target = [gamma{1}(theta,rho) gamma{2}(theta,rho)];
target = reshape(target,[],2);

[M,~] = size(target); % num of target


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Discrete in Gamma (A*mu=f)
% 1.1 b
f = u_exact(src(:,1),src(:,2));

% 1.2 A
eta = [gamma_dt{2}(theta,1) -gamma_dt{1}(theta,1) ];

A = [src(:,1)-src(:,1)' src(:,2)-src(:,2)'];
A = (A(:,1:N).*eta(:,1)' + A(:,N+1:end).*eta(:,2)')./(A(:,1:N).^2+A(:,N+1:end).^2);
A(1:N+1:end) = 1/2*(gamma_dt2{1}(theta).*eta(:,1)+gamma_dt2{2}(theta).*eta(:,2))./(sum(eta.^2,2));
A = 1/N*A-1/2*eye(N);

% 1.3 mu
mu = A\f;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Points in Omega by direct FMM
iprec = 4;
nsource = N;
source = src';
ifcharge = 0;
charge = [];
ifdipole = 1;
dipstr = -1/N*mu;
dipvec = eta';
ifpot = 0;
ifgrad = 0;
ifhess = 0;
ntarget = M;
targets = target';
ifpottarg = 1;
ifgradtarg = 0;
ifhesstarg = 0;

[U] = rfmm2dpart(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,...
        ntarget,targets,ifpottarg,ifgradtarg,ifhesstarg);
if U.ier ~= 0
    error('Error in FMM')
end
u = U.pottarg'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Plot
gamma_plot = @(x,y) x.^2/a^2+y.^2/b^2;
axis_plot = [-a a;-b b];
N_plot = [500,500];
colorbar_range = -15:2:0;
N_contourf = length(colorbar_range);
[u_plot] = myplot(target, log10(abs(u-u_exact(target(:,1),...
    target(:,2)))./u_exact(target(:,1),target(:,2))+1e-15), ...
    gamma_plot, axis_plot,N_plot,N_contourf,colorbar_range);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Points in Omega by correction
% 4.1 partition
rho_reduced = (1-6*h)^2;
target0 = (gamma_plot(target(:,1),target(:,2))<=rho_reduced);
target1 = [target(~target0,1) target(~target0,2)];
target1 = reshape(target1,[],2); % near the boundary
target0 = [target(target0,1) target(target0,2) ];
target0 = reshape(target0,[],2); % away from boundary
[M_target0,~] = size(target0);

% 4.2 target0
ntarget0 = M_target0;
target0s = target0';

[U] = rfmm2dpart(iprec,nsource,source,ifcharge,charge,...
        ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,...
        ntarget0,target0s,ifpottarg,ifgradtarg,ifhesstarg);
if U.ier ~= 0
    error('Error in FMM')
end
u0_correction = U.pottarg'; 


% 4.4 target1
% 4.4.1 paremeters
beta = 4; % interpolation multiple
N_Taylor = floor(beta*N)+1;% interpolation num
theta_Taylor = ((0:N_Taylor-1)/N_Taylor*2*pi)'; % interpolation vector
h_Taylor = 2*pi/N_Taylor; % interpolation h

N_box = floor(N/2)+1; % num of box
p = 15; % Taylor num

theta_z0 = (0:N_box-1)'/N_box*2*pi; % z0
z0 = gamma{1}(theta_z0,1-4*h)+sqrt(-1)*gamma{2}(theta_z0,1-4*h);
% z0 = [z0; gamma{1}(theta_z0,1-3*h)+sqrt(-1)*gamma{2}(theta_z0,1-3*h)];
cm = zeros(p,length(z0)); % coefficient matrix

% 4.4.2 Taylor coefficient
% interpolation
% theta_spline = [(-10:1:-1)'*h;theta;2*pi + (0:9)'*h];
% mu_spline    = [mu(end-9:1:end);mu;mu(1:10)];
% mu_spline    =  spline(theta_spline,mu_spline,theta_Taylor);
mu_spline      =  interpft(mu,N_Taylor);

% zn and tangential
zn = gamma{1}(theta_Taylor,1)+sqrt(-1)*gamma{2}(theta_Taylor,1);
nu = gamma_dt{1}(theta_Taylor,1)+sqrt(-1)*gamma_dt{2}(theta_Taylor,1);
nu = nu.*mu_spline;

% compute
for i=1:p
    cm0 = (zn-z0.').^(-i).*nu;
    cm(i,:) = sum(cm0)*sqrt(-1)/N_Taylor;
end
cm = cm.';

% 4.4.3 value
target1_complex = target1(:,1)+sqrt(-1)*target1(:,2);
[~, index] = min(abs(target1_complex-z0.'), [],2);
target1_2_z0 = target1_complex - z0(index);

u1_correction = cm(index,end); 
for i=p-1:-1:1
    u1_correction = u1_correction.*target1_2_z0 + cm(index,i);
end
u1_correction = real(u1_correction);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Plot
u_correction = [u0_correction;real(u1_correction)];
target_correction = [target0;target1];
[u_correction_plot] = myplot(target, log10(abs(u_correction-u_exact(target_correction(:,1),...
    target_correction(:,2)))./u_exact(target_correction(:,1),target_correction(:,2))+1e-15), ...
    gamma_plot, axis_plot,N_plot, N_contourf,colorbar_range);
hold on
scatter(real(z0),imag(z0),'r')



