%{
This is an implementation of the model described by Dupuy and Fivel

Implementation for a Lomer Lock
%}


clc; clear;

% Definition of important anonymous functions
Energy = @(theta, nu) (1 - nu .* cos(theta).^2) ./ (1-nu);

dEnergydTheta = @(theta, nu) (nu .* sin(2 .* theta)) ./ (1 - nu);

d2Ed2Theta = @(theta, nu) (2.*nu.*cos(2.*theta)) ./ (1 - nu);

Force = @(Psi, Phi, nu) E(Psi - Phi, nu) .* cos(Psi - 2.*Phi) - ...
    dEnergydTheta(Psi - Phi) .* sin(Psi - 2.*Phi)


f = @(Psi_c, Phi_1, Phi_2, Phi_3, nu) Energy(Phi_1 + Psi_c, nu) .* cos(Psi_c) - ...
    dEnergydTheta(Phi_1 + Psi_c, nu) .* sin(Psi_c) + ...
    ...
    Energy(Phi_2 + Psi_c, nu) .* cos(Psi_c) - ...
    dEnergydTheta(Phi_2 + Psi_c, nu) .* sin(Psi_c) - ...
    ...
    Energy(Phi_3, nu);



clc; clear;

% Definition of important anonymous functions
Energy = @(theta, nu) (1 - nu .* cos(theta).^2) ./ (1-nu);
dEnergydTheta = @(theta, nu) (nu .* sin(2 .* theta)) ./ (1 - nu);
d2Ed2Theta = @(theta, nu) (2.*nu.*cos(2.*theta)) ./ (1 - nu);

f = @(Psi_c, Phi_1, Phi_2, Phi_3, nu) Energy(Phi_1 + Psi_c, nu) * cos(Psi_c) + ...
    dEnergydTheta(Phi_1 + Psi_c, nu) * sin(Psi_c) + ...
    ...
    Energy(Phi_2 + Psi_c, nu) * cos(Psi_c) + ...
    dEnergydTheta(Phi_2 + Psi_c, nu) * sin(Psi_c) - ...
    ...
    Energy(Phi_3, nu);

dfdTheta = @(Psi_c, Phi_1, Phi_2, Phi_3, nu) ...
    dEnergydTheta(Phi_1 + Psi_c, nu) * cos(Psi_c) - ...
    Energy(Phi_1 + Psi_c, nu) * sin(Psi_c) + ...
    d2Ed2Theta(Phi_1 + Psi_c, nu) * sin(Psi_c) + ...
    dEnergydTheta(Phi_1 + Psi_c, nu) * cos(Psi_c) + ...
    ...
    dEnergydTheta(Phi_2 + Psi_c, nu) * cos(Psi_c) - ...
    Energy(Phi_2 + Psi_c, nu) * sin(Psi_c) + ...
    d2Ed2Theta(Phi_2 + Psi_c, nu) * sin(Psi_c) + ...
    dEnergydTheta(Phi_2 + Psi_c, nu) * cos(Psi_c) - ...
    ...
    dEnergydTheta(Phi_3, nu);
    
L3L = @(theta, critical) sind(theta ./ abs(theta) .* critical - ...
    theta) ./ sind(theta ./ abs(theta) .* critical);


% Definition of Parameters:
save_loc = './LomerLock/';
nu = 0.3354;

N1 = [1; 1; 1] / norm([1;1;1]);
N2 = [-1; -1; 1] / norm([-1; -1; 1]);

B1 = [0; -1; 1] / norm([0; -1; 1]);
B2 = [-1; 1; 0] / norm([-1; 1; 0]);
B3 = [-1; 0; 1] / norm([-1; 0; 1]);

T_junc = cross(N2, N1);
T_junc = T_junc / norm(T_junc);

Phi1 = -1 * acos(dot(B1, T_junc)); % These are the angles between the burgers vector and the junction direction
Phi2 = acos(dot(B2, T_junc));
Phi3 = -1 * acos(dot(B3, T_junc));

% Calculating the Critical Angle
psi_crit = zeros(1,2);

for i=1:2
    Psi_crit_0 = ((-1)^(2-i)) * 100 * pi / 280;
    psi_crit(i) = newton_raphson(Psi_crit_0, f, dfdTheta, 1e-3, Phi1, Phi2, Phi3, nu);
end
% Plotting the Junction Length as a function of angle using the derived
% formulas

psi_crit
save(append(save_loc, 'psi_criticals.mat'), 'psi_crit');


function final = newton_raphson(initial, f, fprime, error, p1, p2, p3, n)

shift = f(initial, p1, p2, p3, n) / fprime(initial, p1, p2, p3, n);
final = initial;

while abs(shift) > error
    final = final - shift;
    shift = f(final, p1, p2, p3, n) / fprime(final, p1, p2, p3, n);
end

final = final * 180 / pi;

end

