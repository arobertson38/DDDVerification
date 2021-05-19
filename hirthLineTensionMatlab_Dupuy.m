%{
This is an implementation of the model described by Dupuy and Fivel

Implementation for a Hirth Lock
%}

clc; clear;

% Definition of important anonymous functions
Energy = @(theta, nu) (1 - nu .* cos(theta).^2) ./ (1-nu);
dEnergydTheta = @(theta, nu) (nu .* sin(2 .* theta)) ./ (1 - nu);
d2Ed2Theta = @(theta, nu) (2.*nu.*cos(2.*theta)) ./ (1 - nu);

f = @(Psi_c, Phi_1, Phi_2, Phi_3, nu, b_1, b_2, b_3) ...
    (b_1.^2) .* (Energy(Phi_1 + Psi_c, nu) .* cos(Psi_c) + ...
    dEnergydTheta(Phi_1 + Psi_c, nu) .* sin(Psi_c)) + ...
    ...
    (b_2 .^ 2) .* (Energy(Phi_2 + Psi_c, nu) .* cos(Psi_c) + ...
    dEnergydTheta(Phi_2 + Psi_c, nu) .* sin(Psi_c)) - ...
    ...
    (b_3 .^ 2) .* Energy(Phi_3, nu);

dfdTheta = @(Psi_c, Phi_1, Phi_2, Phi_3, nu, b_1, b_2, b_3) ...
    (b_1.^2) .* (dEnergydTheta(Phi_1 + Psi_c, nu) * cos(Psi_c) - ...
    Energy(Phi_1 + Psi_c, nu) * sin(Psi_c) + ...
    d2Ed2Theta(Phi_1 + Psi_c, nu) * sin(Psi_c) + ...
    dEnergydTheta(Phi_1 + Psi_c, nu) * cos(Psi_c)) + ...
    ...
    (b_2 .^ 2) .* (dEnergydTheta(Phi_2 + Psi_c, nu) * cos(Psi_c) - ...
    Energy(Phi_2 + Psi_c, nu) * sin(Psi_c) + ...
    d2Ed2Theta(Phi_2 + Psi_c, nu) * sin(Psi_c) + ...
    dEnergydTheta(Phi_2 + Psi_c, nu) * cos(Psi_c)) - ...
    ...
    (b_3 .^ 2) .* dEnergydTheta(Phi_3, nu);
    
L3L = @(theta, critical) sind(theta ./ abs(theta) .* critical - ...
    theta) ./ sind(theta ./ abs(theta) .* critical);


% Definition of Parameters:
save_loc = './HirthLock/';
nu = 0.347;

N1 = [1; 1; 1] / norm([1;1;1]);
N2 = [1; -1; 1] / norm([1; -1; 1]);

B1 = 0.5 * [0; -1; 1];
B2 = 0.5 * [0;  1; 1];
B3 = B1 + B2;

b1 = norm(B1);
b2 = norm(B2);
b3 = norm(B3);

B1 = B1 / b1;
B2 = B2 / b2;
B3 = B3 / b3;

T_junc = cross(N1, N2);
T_junc = T_junc / norm(T_junc);

Phi1 = acos(dot(B1, T_junc)); % These are the angles between the burgers vector and the junction direction
Phi2 = acos(dot(B2, T_junc));
Phi3 = acos(dot(B3, T_junc));

% Calculating the Critical Angle
psi_crit = zeros(1,2);
pc0 = [-30, 10];

for i=1:2
    Psi_crit_0 = pc0(i) * pi / 180;
    psi_crit(i) = newton_raphson(Psi_crit_0, f, dfdTheta, 1e-3, Phi1, Phi2, Phi3, nu, b1, b2, b3);
end
% Plotting the Junction Length as a function of angle using the derived
% formulas

psi_crit
%save(append(save_loc, 'psi_criticals.mat'), 'psi_crit');

figure(1)
title('Force Function')
x = linspace(-pi/2, pi/2, 1000);
y = f(x, Phi1, Phi2, Phi3, nu, b1, b2, b3);
plot(x,y)


function final = newton_raphson(initial, f, fprime, error, p1, p2, p3, n, b_1, b_2, b_3)

shift = f(initial, p1, p2, p3, n, b_1, b_2, b_3) / fprime(initial, p1, p2, p3, n, b_1, b_2, b_3);
final = initial;

while abs(shift) > error
    final = final - shift;
    shift = f(final, p1, p2, p3, n, b_1, b_2, b_3) / fprime(final, p1, p2, p3, n, b_1, b_2, b_3);
end

final = final * 180 / pi;

end

