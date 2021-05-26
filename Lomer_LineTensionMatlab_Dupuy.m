%{
This is an implementation of the line tension model for predicting the critical
angles to form symmetric dislocation junctions described by Dupuy and Fivel

Implementation for a Lomer Lock
%}

clc; clear;

% Definition of important anonymous functions
Energy = @(theta, b, nu) (b.^2) .* (1 - nu .* cos(theta).^2) ./ (1-nu);

dEnergydTheta = @(theta, b, nu) (b.^2) .* (nu .* sin(2 .* theta)) ./ (1 - nu);

d2Ed2Theta = @(theta, b, nu) (b.^2) .* (2.*nu.*cos(2.*theta)) ./ (1 - nu);

force = @(Psi_c, Phi, b, nu) Energy(Psi_c - Phi, b, nu) .* cos(Psi_c) - ...
    dEnergydTheta(Psi_c - Phi, b, nu) .* sin(Psi_c);

dforcedtheta = @(Psi_c, Phi, b, nu) dEnergydTheta(Psi_c - Phi, b, nu) .* cos(Psi_c) - ...
    Energy(Psi_c - Phi, b, nu) .* sin(Psi_c) - ...
    d2Ed2Theta(Psi_c - Phi, b, nu) .* sin(Psi_c) - ...
    dEnergydTheta(Psi_c - Phi, b, nu) .* cos(Psi_c);



% Combined Functions

f = @(Psi_c, Phi_1, Phi_2, Phi_3, b_1, b_2, b_3, nu) force(Psi_c, Phi_1, b_1, nu) + ...
    force(Psi_c, Phi_2, b_2, nu) - ...
    force(0.0, Phi_3, b_3, nu);


dfdTheta = @(Psi_c, Phi_1, Phi_2, Phi_3, b_1, b_2, b_3, nu) ...
    dforcedtheta(Psi_c, Phi_1, b_1, nu) + ...
    dforcedtheta(Psi_c, Phi_2, b_2, nu) - ...
    dforcedtheta(0.0, Phi_3, b_3, nu);


% -------------------------------------------------------------------


% Defining the Parameters
save_loc = './LomerLock/';
nu = 0.317;

N1 = [1; 1; 1] / norm([1;1;1]);
N2 = [-1; -1; 1] / norm([-1; -1; 1]);

B1 = 0.5 * [0; -1; 1];
B2 = 0.5 * [-1; 1; 0];
B3 = B1 + B2;

b1 = norm(B1);
b2 = norm(B2);
b3 = norm(B3);

B1 = B1 / b1;
B2 = B2 / b2;
B3 = B3 / b3;

T_junc = cross(N1, N2);
T_junc = -1 * T_junc / norm(T_junc);

Phi1 = acos(dot(B1, T_junc)); % These are the angles between the burgers vector and the junction direction
Phi2 = acos(dot(B2, T_junc));
Phi3 = acos(dot(B3, T_junc));


% -------------------------------------------------------------------


% Calculating the Critical Angle
psi_crit = zeros(1,2);

for i=1:2
    Psi_crit_0 = ((-1)^(2-i)) * 60 * pi / 180;
    psi_crit(i) = newton_raphson(Psi_crit_0, f, dfdTheta, 1e-3, Phi1, Phi2, Phi3, b1, b2, b3, nu);
end


% -------------------------------------------------------------------


% Plotting the Junction Length as a function of angle using the derived

%save(append(save_loc, 'psi_criticals.mat'), 'psi_crit');

psi_crit

figure()
title('Force on Triple Point along Junction Direction')
x = linspace(-pi/2, pi/2, 1000);
y = f(x, Phi1, Phi2, Phi3, b1, b2, b3, nu);
plot(x*180/pi,y)
ylabel('Force')
xlabel('\Psi (deg)')
hold on
plot(linspace(-100, 100, 1000), zeros(1,1000), 'k-')

function final = newton_raphson(initial, f, fprime, error, p1, p2, p3, b_1, b_2, b_3, n)

shift = f(initial, p1, p2, p3, b_1, b_2, b_3, n) / fprime(initial, p1, p2, p3, b_1, b_2, b_3, n);
final = initial;

while abs(shift) > error
    final = final - shift;
    shift = f(final, p1, p2, p3, b_1, b_2, b_3, n) / fprime(final, p1, p2, p3, b_1, b_2, b_3, n);
end

final = final * 180 / pi;

end

