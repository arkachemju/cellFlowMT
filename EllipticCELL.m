% Parameters
D = 0.5;  % Diffusion coefficient
dt = 0.01; % Time step
T = 1;    % Total simulation time
a = 6;    % Semi-major axis of the ellipse
b = 4;    % Semi-minor axis of the ellipse

% Define grid
Nx = 100; % Number of points along x-axis
Ny = 100; % Number of points along y-axis
x = linspace(-a, a, Nx);
y = linspace(-b, b, Ny);
[X, Y] = meshgrid(x, y);

% Elliptical boundary condition
ellipse_mask = (X.^2 / a^2 + Y.^2 / b^2) <= 1;

% Initialize Z within the elliptical boundary
Z = zeros(Nx, Ny);
Z(ellipse_mask) = sin(sqrt(X(ellipse_mask).^2 + Y(ellipse_mask).^2));

% Plot the initial mesh
figure;
mesh(X, Y, Z);
title('Initial Meshed Results within Elliptical Boundary');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
shading interp

% Solving the Fokker-Planck equation
for t = 0:dt:T
    % Compute the Laplacian (finite difference approximation)
    Laplacian = (circshift(Z, [0, 1]) + circshift(Z, [0, -1]) + ...
                 circshift(Z, [1, 0]) + circshift(Z, [-1, 0]) - 4 * Z) / (a * b);

    % Update the solution using the Fokker-Planck equation
    Z_new = Z + D * Laplacian * dt;

    % Apply elliptical boundary conditions
    Z = Z_new;
    Z(~ellipse_mask) = 0;
end

% Plotting the results within an elliptical boundary
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
title('Surfed Results within Elliptical Boundary');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis equal;
shading interp

% Overlay ellipse outline
hold on;
theta = linspace(0, 2 * pi, 100);
x_ellipse = a * cos(theta);
y_ellipse = b * sin(theta);
plot3(x_ellipse, y_ellipse, max(Z(:)) * ones(size(theta)), 'r', 'LineWidth', 2);