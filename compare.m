clc
clear
close all
addpath(genpath('..'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 20:5:100;
m= 400;

k = 4;
a = 1+1/(k+1);
b = 1+1/(k+1);
u_exact = @(x,y) 10+1/(a*b)*x.*y +...
                 1/exp(b)*exp(-y).*sin(x) + 1/exp(a)*exp(x).*cos(y);

f = @(location, state) u_exact(location.x,location.y);

cost_time4 = zeros(length(N),2);
error4 = zeros(length(N),2);


for i=1:length(N)
    [target,u_correction,cost_time4(i,1)] = polar_function(N(i),m,k);
    error4(i,1) = norm(abs(u_exact(target(:,1),target(:,2))-u_correction...
        )./u_exact(target(:,1),target(:,2)))/(length(target))*2;
    
    tic
    model = createpde();
    geometryFromEdges(model,@polar4); 
    applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', f);
    specifyCoefficients(model, 'm', 0,'d', 0, 'c', 1, 'a', 0, 'f', 0);
    mesh_default = generateMesh(model,'Hmax',1/N(i));
    result = solvepde(model);
    cost_time4(i,2) = toc;
    result4 = interpolateSolution(result, target(:,1), target(:,2));
    error4(i,2) = norm(abs(u_exact(target(:,1),target(:,2))-result4...
        )./u_exact(target(:,1),target(:,2)))/(length(target))*2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 5;
a = 1+1/(k+1);
b = 1+1/(k+1);

u_exact = @(x,y) 10+1/(a*b)*x.*y +...
                 1/exp(b)*exp(-y).*sin(x) + 1/exp(a)*exp(x).*cos(y);

f = @(location, state) u_exact(location.x,location.y);

cost_time5 = zeros(length(N),2);
error5 = zeros(length(N),2);

for i=1:length(N)
    [target,u_correction,cost_time5(i,1)] = polar_function(N(i),m,k);
    error5(i,1) = norm(abs(u_exact(target(:,1),target(:,2))-u_correction...
        )./u_exact(target(:,1),target(:,2)))/(length(target))*2;
    
    tic
    model = createpde();
    geometryFromEdges(model,@polar5); 
    applyBoundaryCondition(model, 'dirichlet', 'Edge', 1:model.Geometry.NumEdges, 'u', f);
    specifyCoefficients(model, 'm', 0,'d', 0, 'c', 1, 'a', 0, 'f', 0);
    mesh_default = generateMesh(model,'Hmax',1/N(i));
    result = solvepde(model);
    cost_time5(i,2) = toc;
    result5 = interpolateSolution(result, target(:,1), target(:,2));
    error5(i,2) = norm(abs(u_exact(target(:,1),target(:,2))-result5...
        )./u_exact(target(:,1),target(:,2)))/(length(target))*2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(N,log10(error4(:,1)),'LineWidth',2)
hold on
plot(N,log10(error4(:,2)),'LineWidth',2)
xlabel('N')
ylabel('log_{10}(error)')
legend('Integral-Method','FEM')

figure(2)
plot(N,log10(cost_time4(:,1)),'LineWidth',2)
hold on
plot(N,log10(cost_time4(:,2)),'LineWidth',2)
xlabel('N')
ylabel('log_{10}(time)')
legend('Integral-Method','FEM')

figure(11)
plot(N,log10(error5(:,1)),'LineWidth',2)
hold on
plot(N,log10(error5(:,2)),'LineWidth',2)
xlabel('N')
ylabel('log_{10}(error)')
legend('Integral-Method','FEM')

figure(21)
plot(N,log10(cost_time5(:,1)),'LineWidth',2)
hold on
plot(N,log10(cost_time5(:,2)),'LineWidth',2)
xlabel('N')
ylabel('log_{10}(time)')
legend('Integral-Method','FEM')




