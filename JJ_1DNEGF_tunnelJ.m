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

N_p = 50;
phi_mat = linspace(0, 2*pi, N_p); 
Ij_mat = zeros(1, N_p);

for k = 1:N_p
    phi = phi_mat(k);k
    ep = 0.02;                                                                   % onsite energy

    alpha_1 = [ep - mu1, del_1*exp(-1i * phi); del_1*exp(1i * phi), -(ep - mu1)];                        % onsite interaction for lead 1
    alpha_2 = [ep - mu2, del_1; del_1, -(ep - mu2)];                        % onsite interaction for lead 2

    beta = [-t, 0; 0, t];                                                       % hopping parameters

    M =0.05 * eye(2);                                                          % Tunneling rate

    I_J = 0;

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

        I_J = I_J + 4 * real(gs1(2,1)*gs2(1,2) - gs2(2,1)*gs1(1,2)) * f1(j);

    end
    Ij_mat(k) = I_J;
end
%% Plots
figure(6)                                                                                                   % Subpart c                                             
h = plot(phi_mat/(2*pi),Ij_mat);
hold on
g = plot(phi_mat/(2*pi), max(Ij_mat) * -sin(phi_mat));
set(h,'linewidth',2.0)
set(gca,'Fontsize',24)
ylabel("current(arb units)");
xlabel("\Phi (rad)");
xlim padded
grid on

%% Quasi-particle current
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar; B=0; kT = 3 * 8.6173303 * 1e-5;       % Lattice Constants
%inputs
a=2.5e-9; t=0.1; N = 2;

zplus=1i*1e-12;                                                        % Infinitesimal eta

% Levels and s-wave junction parameters
mu_f = 0; del_d = t * 0.05; del0 = 0.0015;                       % Pairing parameters
del_1 = del_d ; del_2 =  del_d * 0.5 ;
Emax = 5 * del0; Emin = - 5 * del0; Ne = 2000; E = linspace(Emin, Emax, Ne);
n_vol = 40;
vol_mat = linspace(-4*(del_1 + del_2), 4*(del_1 + del_2), n_vol);

IQP_mat = zeros(1, n_vol);
Ij_mat2 =  zeros(1, n_vol);


ppm = ParforProgressbar(n_vol);
parfor k = 1:n_vol
    vol = vol_mat(k);k
%     mu1 = mu_f + vol/2; mu2 = mu_f - vol/2;

    f1=1./(1+exp((E + vol/2  - mu_f)./kT));
    f2=1./(1+exp((E  - vol/2 - mu_f)./kT));
    ep = 0.02;                                                                   % onsite energy
    phi = 0;
    alpha_1 = [ep - mu_f, del_1*exp(-1i * phi); del_1*exp(1i * phi), -(ep - mu_f)];                        % onsite interaction for lead 1
    alpha_2 = [ep - mu_f, del_2; del_2, -(ep - mu_f)];                        % onsite interaction for lead 2

    beta = [-t, 0; 0, t];                                                       % hopping parameters

    M =0.05 * eye(2);                                                          % Tunneling rate

    I_J = 0;IQP = 0;

    for j = 1:Ne
        gs1=inv((E(j) + vol/2  + zplus)*eye(2,2));gs2=inv((E(j) - vol/2 + zplus)*eye(2,2));
        % Green's function for Lead 1
        change=1;
        while change >1e-4
            g1=inv((E(j) + vol/2 + zplus)*eye(2,2) - alpha_1  -beta'*gs1*beta);
            change=sum(sum(abs(g1-gs1)))/(sum(sum(abs(g1)+abs(gs1))));
            gs1=0.5*g1+0.5*gs1;
        end
        % Green's function for Lead 2
        change=1;
        while change >1e-4
            g2=inv((E(j) -vol/2 + zplus)*eye(2,2) - alpha_2  -beta*gs2*beta');
            change=sum(sum(abs(g2-gs2)))/(sum(sum(abs(g2)+abs(gs2))));
            gs2=0.5*g2+0.5*gs2;
        end

        I_J = I_J + 4 * real(gs1(2,1)*gs2(1,2) - gs2(2,1)*gs1(1,2)) * f1(j);

        a_1 = 1i * (gs1 - gs1'); a_2 = 1i * (gs2 - gs2');

        IQP = IQP + (a_1(1,1)*a_2(1,1) + a_1(2,2)*a_2(2,2)) * (f2(j) - f1(j));
    end
    IQP_mat(k) = IQP;Ij_mat2(k) = I_J;
    pause(100/n_vol);
    ppm.increment(); 
end
delete(ppm)
%% Plots
figure(7)                                                                                                
hold on
h = plot(vol_mat/(del_1 + del_2),IQP_mat, 'DisplayName',...
    "\Delta_1 = " + del_1 + ", \Delta_2 = " + del_2);
set(h,'linewidth',2.0)
set(gca,'Fontsize',24)
ylabel("current(arb units)");
xlabel("eV/ (\Delta_1 + \Delta_2)");
legend(gca,'show')
xlim padded
grid on


