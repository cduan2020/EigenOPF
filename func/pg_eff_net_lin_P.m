function [P, x0] = pg_eff_net_lin_P(ps, YEN)
% pg_eff_net_lin_P calculates the P matrix from a given ps and effective
% admittances. Optionally, the linearization point is also returned.

%% Get constants
x_d = ps.gen_dyn(:,1); %transient reactances
H = ps.gen_dyn(:,2);   %inertia constants
baseMVA = ps.baseMVA;  %common system base

% System reference frequency omega_R (in radian):
if isfield(ps, 'ref_freq')
    omega_R = 2 * pi * ps.ref_freq;
else
    omega_R = 2 * pi * 60;    
end

%% calculate voltages at generator internal nodes

% P and Q injected at generator terminals in p.u. on system base MVA.
Pi = ps.gen(:,2) / baseMVA; 
Qi = ps.gen(:,3) / baseMVA;

% Voltage magnitude V and phase angle phi for the generator terminal buses
tb =  ps.gen(:, 1);
V =   ps.bus(tb, 8);
phi = ps.bus(tb, 9) / 180 * pi;

% Compute the complex voltage E at the internal nodes of generators and
% motors.
E = ((V + Qi .* x_d ./ V) + 1j * (Pi .* x_d ./ V)) .* exp(1j * phi);

%% Calculate P

delta = angle(E);
E_abs = abs(E);
%G = real(YEN);
%B = imag(YEN);

Yabs = abs(YEN);
gamma = angle(YEN);

n = size(YEN,1);

% faster this way - check needed
DIK = bsxfun(@minus, delta, delta.'); %delta(i) - delta(k)
K = bsxfun(@times, E_abs, E_abs') .* Yabs; % |E(i)*E(k)*Y(i,k)|
P = -K .* cos(DIK - gamma + pi/2); 
% for i = 1 : n
%     P(i,i) = - sum(P(i,[1:i-1,i+1:n]));
%     P(i,:) = P(i,:) * omega_R / (2*H(i));
% end

P0=P;
PP=P0-diag(diag(P0))-diag(sum(P0-diag(diag(P0)),2));
PP=PP.*((omega_R./(2*H)).*ones(1,n));
P=PP;

P=diag(sqrt(H))*P*diag(1./sqrt(H));
% norm(full(P-PP))

% tic;
% % Compute P
% ddik = zeros(n,n);
% for i = 1:n
%     for k = [1:i-1, i+1:n]
%         dik = delta(i) - delta(k);
%         ddik(i,k) = dik;
%         P(i,k) = - E_abs(i)*E_abs(k)/(2*H(i))*(B(i,k)*cos(dik) - G(i,k)*sin(dik));
%     end
% end
% for i = 1:n
%     P(i,i) = - sum(P(i,[1:i-1,i+1:n]));
% end
% P = omega_R * P;
% % -------check here: P2 == P -- yes. normdiff=1e-13
% toc;

x0 = delta;

end