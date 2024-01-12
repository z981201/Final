function [u_plot] = myplot(target, u, gamma_plot, axis_plot,N_plot,N_contourf,colorbar_range)
% Input:
% target: target's coordinate, (M,2)-vector
% u    :  target's value     ,(M,1)-vector
% gamma_plot:  decide interior point    , =f(x,y)(<1)
% axis_plot: plot range      ,[x_min, x_max; y_min, y_max]
% N_plot:  plot mesh         ,[x_N, y_N]
% Outpu:
% u_plot: figure            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot mesh
[X, Y] = meshgrid(linspace(axis_plot(1,1),axis_plot(1,2),N_plot(1)), ...
                      linspace(axis_plot(2,1),axis_plot(2,2),N_plot(2)));
Z = griddata(target(:,1), target(:,2), u, X, Y, 'cubic');

% interior points and boundary points
interior_points = gamma_plot(X,Y)<=1;
Z(~interior_points) = nan;

% plot
u_plot = figure;
contourf(X, Y, Z, N_contourf);  
hold on;
contour(X, Y, Z, N_contourf,'k');  
% colorbar;
caxis([colorbar_range(1) colorbar_range(end)]);
cb = colorbar;
cb.Ticks = colorbar_range;
axis equal;
axis off;
hold off;
end
