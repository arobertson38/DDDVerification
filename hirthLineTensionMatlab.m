%{
Trying to implement a line tension model in Matlab


%}

% Defining necessary variables

Energy = @(theta, nu) (1 - nu .* cos(theta).^2) ./ (1-nu);
dEnergydTheta = @(theta, nu) (nu .* sin(2 .* theta)) ./ (1 - nu);

DeltaEnergy = @(theta1, theta2, theta3, b1, b2, b3, phi1, phi2, nu) ...
    (b3^2 * (1-nu*cos(theta3)^2)) - (b1^2 * (1-nu*cos(theta1)^2) * cos(phi1) ...
    ) - (b2^2 * (1-nu * cos(theta2)^2) * cos(phi2));

b_1 = [0; -1; 1] / norm([0; -1; 1]);
n_1 = [1; 1; 1] / norm([1;1;1]);
e_1 = cross(b_1, n_1);

b_2 = [0; 1; 1] / norm([0; 1; 1]);
n_2 = [1; -1; 1] / norm([1; -1; 1]);
e_2 = cross(b_2, n_2);


b_3 = [0; 0; 1] / norm([0; 0; 1]);
t_3 = cross(n_1, n_2);
theta_3 = acos(dot(b_3, t_3));

L = 400;
nu = 0.3; 

b1_mag = (0.5^2 + 0.5^2)^(1/2);
b2_mag = (0.5^2 + 0.5^2)^(1/2);
b3_mag = 1;

% calculation:
Lens = [];
dEnergies = zeros(20,20); %[];

% Extra
count_1 = 1;
count_2 = 1;
%

for theta_1 = linspace(0, 2*pi, 20)
    for theta_2 = linspace(0, 2*pi, 20)
        % Calculating the Length
        t_1 = sin(theta_1) * e_1 + cos(theta_1) * b_1;
        t_2 = sin(theta_2) * e_2 + cos(theta_2) * b_2;
        
        F_1A = (L/2) * (Energy(theta_1, 0.3) * t_1 + dEnergydTheta(theta_1, ...
            0.3) * n_1);
        F_2A = (L/2) * (Energy(theta_2, 0.3) * t_2 + dEnergydTheta(theta_2, ...
            0.3) * n_2);
        
        L_3 = dot(F_1A + F_2A, t_3) ./ Energy(theta_3, 0.3);
        Lens = [Lens, L_3 / L];
        
        % Calculating the energy change due to combination
        phi1 = acos(dot(t_3, t_1));
        phi2 = acos(dot(t_3, t_2));
        
        dEnergies(count_2, count_1) = DeltaEnergy(theta_1, theta_2, theta_3, ...
            b1_mag, b2_mag, b3_mag, phi1, phi2, nu); %[dEnergies, DeltaEnergy(theta_1, theta_2, theta_3, ...
            %b1_mag, b2_mag, b3_mag, phi1, phi2, nu)];
        count_2 = count_2 + 1;
    end
    count_1 = count_1 + 1;
    count_2 = 1;
end

[x,y] = meshgrid(linspace(0, 2*pi, 20), linspace(0, 2*pi, 20));
figure();
contour(x,y,dEnergies, 'ShowText', 'on');
xlabel('\Theta_1');
ylabel('\Theta_2');
title('\Delta E / D')

% Plotting Length
%figure();
%plot(linspace(0, 2*pi, 20) * 180 / pi, Lens, 'ks-');
%xlabel('\Theta_1');
%ylabel('L_3 / L');
%title('\Theta_2 = \pi/2')

% Plotting Energy
%figure(2);
%plot(linspace(0, 2*pi, 20) * 180 / pi, dEnergies, 'kd-')
%xlabel('\Theta_1')
%ylabel('\Delta E / D')
%title('\Theta_2 = \pi / 2')




    
