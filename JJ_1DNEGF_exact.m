clear all
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar; B=0; kT = 0.025;                    % Lattice Constants
%inputs
a=2.5e-9; t=0.1; N = 2;

zplus=1i*1e-12;                                                        % Infinitesimal eta

% Levels and s-wave junction parameters
mu_f = 1.5 * t; del_d = t * 0.05; del0 = del_d * mu_f/t;                       % Pairing parameters
del_1 = del_d * 0.9; del_2 =  del_d * 0.8;
Emax = 5 * del0; Emin = - 5 * del0; Ne = 30; E = linspace(Emin, Emax, Ne);
vol = 0.00; 
mu1 = mu_f + vol/2; mu2 = mu_f + -vol/2;

f1=1./(1+exp((E-mu1)./kT));
f2=1./(1+exp((E-mu2)./kT));

N_p = 30;
phi_mat = linspace(0, 2*pi, N_p); 
Ij_mat_exact = zeros(1, N_p);
Ij_mat_approx = zeros(1, N_p);

for k = 1:N_p
    phi = phi_mat(k);k
    ep = 0.02;                                                                   % onsite energy

    alpha_1 = [ep - mu1, del_1*exp(-1i * phi); del_1*exp(-1i * phi), -(ep - mu1)];                        % onsite interaction for lead 1
    alpha_2 = [ep - mu2, del_1; del_1, -(ep - mu2)];                        % onsite interaction for lead 2

    beta = [-t, 0; 0, t];                                                       % hopping parameters

    M =0.05 * eye(2);                                                          % Tunneling rate

    I_J_ex = 0;
    I_J_app = 0;

    for j = 1:Ne
        gs1=inv((E(j) + zplus)*eye(2,2));gs2=inv((E(j) + zplus)*eye(2,2));
        % Green's function for Lead 1
        change=1;
        while change >1e-6
            g1=inv((E(j) + zplus)*eye(2,2) - alpha_1  -beta'*gs1*beta);
            change=sum(sum(abs(g1-gs1)))/(sum(sum(abs(g1)+abs(gs1))));
            gs1=0.5*g1+0.5*gs1;
        end
        % Green's function for Lead 2
        change=1;
        while change >1e-6
            g2=inv((E(j) + zplus)*eye(2,2) - alpha_2  -beta*gs2*beta');
            change=sum(sum(abs(g2-gs2)))/(sum(sum(abs(g2)+abs(gs2))));
            gs2=0.5*g2+0.5*gs2;
        end

        G = inv([inv(gs1) -M; -M', inv(gs2)]);
        sig_1 = beta' * gs1 * beta; sig_2 = beta * gs2 * beta';
 
        sigLess_1 = 1i * (sig_1 - sig_1') * f1(j); sigLess_2 = 1i * (sig_2 - sig_2') * f2(j);
        G_less =  G * ([sigLess_1, zeros(2,2); zeros(2,2), sigLess_2]) * G';
 
        tau_3 = [1 0; 0 -1];

        % I_op = q/(hbar * 2* pi) * det(M)^2 * (tau_3 * gs1 *tau_3 * alpha_1 + tau_3 *alpha_2 * tau_3 * gs2 - ...
        %     tau_3*alpha_1*tau_3*gs2' - tau_3*gs1*tau_3*alpha_2); 
        I_J_ex = (q/(hbar * 2* pi)) * trace(M * G_less(1:2, 3:4) - G_less(3:4, 1:2) * M);

        I_J_ex = I_J_ex + 4 * real(G(3,1) - G(1,3) + G(4,2) - G(2,4)) * f1(j);

        I_J_app = I_J_app + 4 * real(gs1(2,1)*gs2(1,2) - gs2(2,1)*gs1(1,2)) * f1(j);

    end
    Ij_mat_exact(k) = I_J_ex;
    Ij_mat_approx(k) = I_J_app;
end
%% Plots
figure(6)                                                                                                   % Subpart c                                             
h = plot(phi_mat/(2*pi),Ij_mat);
hold on
g = plot(phi_mat/(2*pi),Ij_mat_approx);
set(h,'linewidth',2.0)
set(g,'linewidth',2.0)
set(gca,'Fontsize',24)
ylabel("current(arb units)");
xlabel("\Phi (rad)");
xlim padded
grid on

