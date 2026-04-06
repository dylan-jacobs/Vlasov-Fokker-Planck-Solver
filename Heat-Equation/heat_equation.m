%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----- Heat Equation Temporal Accuracy Test ----- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear variables; close all;

% initial rank
r0 = 10;

% time-stepping method: 1=B.Euler, 2=DIRK2, 3=DIRK3
methods = ['1', '2', '3', '4'];
tolerance = 1e-6;

lambdavals = (0.02:0.02:10)'; % lambdavals = 0.1;
errors = zeros(numel(lambdavals), 4, 2);

Nrvals = [40, 80];
Nzvals = [80, 160];

for spatial_accuracy_test_idx = 1:1

    % mesh parameters
    rmin = 0; rmax = 1;
    zmin = 0; zmax = 1;
    % Nr = 40; Nz = 80;
    Nr = Nrvals(spatial_accuracy_test_idx);
    Nz = Nzvals(spatial_accuracy_test_idx);
    [Rmat, Zmat, dr, dz] = GetRZ(rmin, rmax, zmin, zmax, Nr, Nz);
    rvals = Rmat(:, 1);
    zvals = Zmat(1, :)';
    
    % final time
    tf = 0.1;
    
    Drr = gallery('tridiag', Nr, rvals(1:end-1)+(dr/2), -2*rvals, rvals(1:end-1)+(dr/2));
    Drr(1, 1) = -(rvals(1) + (dr/2));
    Drr(end, end-1) = ((1/3) * (rvals(end) + (dr/2))) + (rvals(end) - (dr/2));
    Drr(end, end) =  ((-3) * (rvals(end) + (dr/2))) - (rvals(end) - (dr/2));
    Drr = (1/(dr^2)) * (diag(1./rvals) * Drr);
    
    Dzz = (1/dz^2)*gallery('tridiag', Nz, 1, -2, 1); % centered nodes
    
    % Dzz = (2*pi/zmax)^2*toeplitz([-1/(3*(2*dz/zmax)^2)-1/6 ...
    % .5*(-1).^(2:Nz)./sin((2*pi*dz/zmax)*(1:Nz-1)/2).^2]);
 
    for x = 1:3
        method = methods(x);
        
        for k = 1:numel(lambdavals)
            dt = lambdavals(k)/((1/dr) + (1/dz));
            disp(['Method: ', method, ', dt=', num2str(dt, 3), ', ', num2str(k), '/', num2str(numel(lambdavals))]);
            tvals = (0:dt:tf)';
            if tvals(end) ~= tf
                tvals = [tvals; tf];
            end
            Nt = numel(tvals);
            
            % % initial conditions
            j01 = 2.40482555769577; % first root of bessel function
            f0 = @(r, z) (besselj(0, (j01/(rmax-rmin))*r)) .* (sin((2*pi/(zmax-zmin))*z));
            f_exact = @(r, z, t) (exp(-t*(((j01/(rmax-rmin))^2) + (2*pi/(zmax-zmin))^2))*f0(r, z));
            f_exact = f_exact(Rmat, Zmat, tf);
            
            f = f0(Rmat, Zmat);
    
            rhoM = sum(sum(f.*Rmat))*2*pi*dr*dz;
            JzM  = sum(sum(f.*Rmat.*Zmat))*2*pi*dr*dz;
            kappaM   = sum(sum(f.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;
            
            % init bases
            [Vr, S, Vz] = svd2(f, rvals);
            r0 = min(r0, size(Vr, 2));
            Vr = Vr(:, 1:r0); S = S(1:r0, 1:r0); Vz = Vz(:, 1:r0);
            
            % store rank, mass, momentum, energy, l1 decay, etc...
            mass = zeros(Nt, 1);
            Jzvals = zeros(Nt, 1);
            E = zeros(Nt, 1);
            min_vals = zeros(Nt, 1);
            ranks = zeros(Nt, 1);
            
            mass(1) = 2*pi*dr*dz*sum(sum(Rmat .* f));
            Jzvals(1) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
            E(1) = 2*pi*dr*dz*sum(sum(f .* Rmat .* ((Rmat.^2 + Zmat.^2)/2)));
            min_vals(1) = min(min(f));
            ranks(1) = r0;
            
            % time-stepping loop
            for n = 2:Nt
                tval = tvals(n);
                dt = tval - tvals(n-1);
                switch(method)
                    case '1'
                        [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                    case '2'
                        [Vr, S, Vz, rank] = DIRK2Timestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                    case '3'
                        [Vr, S, Vz, rank] = DIRK3Timestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                    case '4'
                        [Vr, S, Vz, rank] = DIRK4Timestep(Vr, S, Vz, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
                end
                
                f = Vr*S*Vz';
                mass(n) = 2*pi*dr*dz*sum(sum(Rmat .* f));    
                Jzvals(n) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
                E(n) = pi*dr*dz*sum(sum(f .* (Rmat.^2 + Zmat.^2) .* Rmat));
                min_vals(n) = min(min(f));
                ranks(n) = rank;
            end
            
            errors(k, x, spatial_accuracy_test_idx) = 2*pi*dr*dz*sum((sum(Rmat .* abs(f - f_exact)))); % L1 error
        end
    end
end

%%
% 1. Final solution
figure(1); clf; surf(Rmat, Zmat, f);
colorbar; shading interp;
legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
xlabel('V_r'); ylabel('V_z'); zlabel('U'); title([sprintf('Backward Euler approximation of 0D2V Fokker-Planck system at time %s', num2str(tf, 4))]);

% 2. Exact solution
figure(2); clf; surf(Rmat, Zmat, f_exact);
colorbar; shading interp;
xlabel('V_r'); ylabel('V_z'); zlabel('f(V_r, V_z, t)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);

%% 3. Temporal error plot
dtvals = lambdavals./((1/dr) + (1/dz));
figure(3); clf; 
begin_cutoff_1 = ceil(0.05*numel(dtvals)); end_cutoff_1 = ceil(0.4*numel(dtvals));
begin_cutoff_2 = ceil(0.2*numel(dtvals)); end_cutoff_2 = ceil(0.7*numel(dtvals));
begin_cutoff_3 = ceil(0.4*numel(dtvals)); end_cutoff_3 = ceil(1*numel(dtvals));

% Nr=40, Nz=80
loglog(dtvals, errors(:, 1, 1), 'black-', 'LineWidth', 1.5); hold on; % B. Euler
loglog(dtvals, errors(:, 2, 1), 'blue-', 'LineWidth', 1.5); % DIRK2
loglog(dtvals, errors(:, 3, 1), 'green-', 'LineWidth', 1.5); % DIRK3

% Nr=80, Nz=160
loglog(dtvals, errors(:, 1, 2), 'black-', 'LineWidth', 1.5); hold on; % B. Euler
loglog(dtvals, errors(:, 2, 2), 'blue-', 'LineWidth', 1.5); % DIRK2
loglog(dtvals, errors(:, 3, 2), 'green-', 'LineWidth', 1.5); % DIRK3

% loglog(dtvals, errors(:, 4), 'm-', 'LineWidth', 1.5); % DIRK4
% loglog(dtvals(begin_cutoff_1:end_cutoff_1), 0.000002*dtvals(begin_cutoff_1:end_cutoff_1), 'black--', 'LineWidth', 1); % Order 1
loglog(dtvals(begin_cutoff_1:end_cutoff_1), 0.8*dtvals(begin_cutoff_1:end_cutoff_1), 'black--', 'LineWidth', 1); % Order 1
% loglog(dtvals(begin_cutoff_2:end_cutoff_2), 0.00000015*dtvals(begin_cutoff_2:end_cutoff_2).^2, 'blue--', 'LineWidth', 1); % Order 2
loglog(dtvals(begin_cutoff_2:end_cutoff_2), 4*dtvals(begin_cutoff_2:end_cutoff_2).^2, 'blue--', 'LineWidth', 1); % Order 2
% loglog(dtvals(begin_cutoff_3:end_cutoff_3), 0.000000008*dtvals(begin_cutoff_3:end_cutoff_3).^3, 'green--', 'LineWidth', 1); % Order 3
loglog(dtvals(begin_cutoff_3:end_cutoff_3), 24*dtvals(begin_cutoff_3:end_cutoff_3).^3, 'green--', 'LineWidth', 1); % Order 3
% loglog(dtvals(ceil(cutoff*end):end), 0.000009*dtvals(ceil(cutoff*end):end).^4, 'm--', 'LineWidth', 1); % Order 4
% title(sprintf('RAIL Temporal Convergence at tf=%s, Nr = %s, Nz = %s', num2str(tf), num2str(Nr), num2str(Nz))); 
xlabel('\Deltat'); ylabel('L_1 Error');
legend('Backward Euler', 'DIRK2', 'DIRK3', '', '', '', 'Order 1', 'Order 2', 'Order 3', 'location','northwest');
fontsize(18,"points");
set(gcf,'Units','pixels','Position',[100 100 800 500])
saveas(gcf, './Plots/heat_eqn_temporal_error_single_spatial.fig');
exportgraphics(gcf,'./Plots/heat_eqn_temporal_error.pdf','ContentType','vector');

%% 4. Mass conservation
figure(4); clf; 
plot(tvals(2:end), abs(mass(2:end)-mass(1))/mass(1), 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative mass'); title('Relative mass of numerical solution');















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- HELPER FUNCTIONS ---- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    Nr = numel(rvals);
    Nz = numel(zvals);

    Vr0_star = Vr0;
    Vz0_star = Vz0;

    K0 = Vr0*S0;
    L0 = Vz0*S0';

    K1 = sylvester(eye(Nr) - (dt*Drr), -dt*(Dzz*Vz0_star)'*Vz0_star, K0);
    L1 = sylvester(eye(Nz) - (dt*Dzz), -dt*(Drr*Vr0_star)'*(rvals.*Vr0_star), L0);

    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);

    S1_hat = sylvester((speye(size(Vr1_hat, 2)) - (dt*((rvals .* Vr1_hat)')*(Drr*Vr1_hat))), -dt*(Dzz*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)'*Vr0)*S0*((Vz0')*Vz1_hat));
    [Vr, S, Vz, rank] = LoMaC_mass_only(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM);
    % [Vr, S, Vz, rank] = truncate_svd(Vr1_hat, S1_hat, Vz1_hat, tolerance);
end

function [Vr, S, Vz, rank] = DIRK2Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    Nr = numel(rvals);
    Nz = numel(zvals);

    nu = 1-(sqrt(2)/2);

    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, nu*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);


    W0 = (Vr0*S0*(Vz0')) + ((1-nu)*dt*(((Drr*(Vr1)*S1*(Vz1')) + ((Vr1)*S1*((Dzz*Vz1)')))));

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vr1_dagger, ~, Vz1_dagger, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_star, Vz_star] = reduced_augmentation([Vr1_dagger, Vr1, Vr0], [Vz1_dagger, Vz1, Vz0], rvals);
   
    % K/L-Step
    K1 = sylvester(eye(Nr) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star)'*Vz_star, W0*Vz_star);
    L1 = sylvester(eye(Nz) - (nu*dt*Dzz), -nu*dt*(Drr*Vr_star)'*(rvals .* Vr_star), W0'*(rvals .* Vr_star));

    % Get bases
    [Vr_ddagger, ~] = qr2(K1, rvals); [Vz_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr_ddagger, Vr1, Vr0], [Vz_ddagger, Vz1, Vz0], rvals);

    % S-Step
    S1_hat = sylvester(eye(size(Vr1_hat, 2)) - (nu*dt*((rvals .* Vr1_hat)')*Drr*Vr1_hat), -nu*dt*(Dzz*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)')*W0*Vz1_hat);
    [Vr, S, Vz, rank] = LoMaC_mass_only(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM);
    % [Vr, S, Vz, rank] = truncate_svd(Vr1_hat, S1_hat, Vz1_hat, tolerance);
end

function [Vr3, S3, Vz3, rank] = DIRK3Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    % RK butcher table values
    nu = 0.435866521508459;
    beta1 = -(3/2)*(nu^2) + (4*nu) - (1/4);
    beta2 = (3/2)*(nu^2) - (5*nu) + (5/4);
    
    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, nu*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_dagger1, ~, Vz_dagger1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, ((1+nu)/2)*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    Y1 = (((Drr*Vr1*S1*(Vz1')) + (Vr1*S1*((Dzz*Vz1)'))));
    W1 = (Vr0*S0*(Vz0')) + (((1-nu)/2)*dt*Y1);

    % Reduced Augmentation
    [Vr_star1, Vz_star1] = reduced_augmentation([Vr_dagger1, Vr1, Vr0], [Vz_dagger1, Vz1, Vz0], rvals);

    % Stage 2:  
    % K/L-Step
    K2 = sylvester(eye(size(Drr)) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star1)'*Vz_star1, W1*Vz_star1);
    L2 = sylvester(eye(size(Dzz)) - (nu*dt*Dzz), -nu*dt*(Drr*(rvals .* Vr_star1))'*Vr_star1, (W1')*(rvals .* Vr_star1));

    % Get bases
    [Vr_ddagger2, ~] = qr2(K2, rvals); [Vz_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vr2, Vz2] = reduced_augmentation([Vr_ddagger2, Vr1, Vr0], [Vz_ddagger2, Vz1, Vz0], rvals);

    % S-Step
    S2 = sylvester(eye(size(Vr2, 2)) - (nu*dt*(rvals .* Vr2)'*Drr*Vr2), -nu*dt*(Dzz*Vz2)'*Vz2, ((rvals .* Vr2)')*W1*Vz2);
    [Vr2, S2, Vz2, ~] = truncate_svd(Vr2, S2, Vz2, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vr_dagger3, ~, Vz_dagger3, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    Y2 = (((Drr*Vr2*S2*(Vz2')) + (Vr2*S2*((Dzz*Vz2)'))));
    W2 = (Vr0*S0*(Vz0')) + (beta1*dt*Y1) + (beta2*dt*Y2);
      
    % Reduced augmentation
    [Vr_star3, Vz_star3] = reduced_augmentation([Vr_dagger3, Vr2, Vr1, Vr0], [Vz_dagger3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K3 = sylvester(eye(size(Drr)) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star3)'*Vz_star3, W2*Vz_star3);
    L3 = sylvester(eye(size(Dzz)) - (nu*dt*Dzz), -nu*dt*(Drr*Vr_star3)'*(rvals .* Vr_star3), (W2')*(rvals .* Vr_star3));

    % Get bases
    [Vr_ddagger3, ~] = qr2(K3, rvals); [Vz_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vr3_hat, Vz3_hat] = reduced_augmentation([Vr_ddagger3, Vr2, Vr1, Vr0], [Vz_ddagger3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S3_hat = sylvester(eye(size(Vr3_hat, 2)) - (nu*dt*((rvals .* Vr3_hat)')*Drr*Vr3_hat), -nu*dt*(Dzz*Vz3_hat)'*Vz3_hat, ((rvals .* Vr3_hat)')*W2*Vz3_hat);
    [Vr3, S3, Vz3, rank] = LoMaC_mass_only(Vr3_hat, S3_hat, Vz3_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM);
    % [Vr3, S3, Vz3, rank] = truncate_svd(Vr3_hat, S3_hat, Vz3_hat, tolerance);
end


function [Vr4, S4, Vz4, rank] = DIRK4Timestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM)

    % RK butcher table values
    a_kk = 0.5;

    % # cvals = np.array([1/2,1/4,3/2,1])
    % 
    % # avals = np.array([[1/2,0,0,0],
    % #                  [-1/4,1/2,0,0],
    % #                  [-1,2,1/2,0],
    % #                  [-1/12,2/3,-1/12,1/2]
    % # ])

    
    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, a_kk*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    [Vr_dagger1, ~, Vz_dagger1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, 0.25*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);

    Y1 = (((Drr*Vr1*S1*(Vz1')) + (Vr1*S1*((Dzz*Vz1)'))));
    W1 = (Vr0*S0*(Vz0')) + (-(1/4)*dt*Y1);

    % Reduced Augmentation
    [Vr_star1, Vz_star1] = reduced_augmentation([Vr_dagger1, Vr1, Vr0], [Vz_dagger1, Vz1, Vz0], rvals);

    % Stage 2:  
    % K/L-Step
    K2 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star1)'*Vz_star1, W1*Vz_star1);
    L2 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*(rvals .* Vr_star1))'*Vr_star1, (W1')*(rvals .* Vr_star1));

    % Get bases
    [Vr_ddagger2, ~] = qr2(K2, rvals); [Vz_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vr2_hat, Vz2_hat] = reduced_augmentation([Vr_ddagger2, Vr1, Vr0], [Vz_ddagger2, Vz1, Vz0], rvals);

    % S-Step
    S2_hat = sylvester(eye(size(Vr2_hat, 2)) - (a_kk*dt*(rvals .* Vr2_hat)'*Drr*Vr2_hat), -a_kk*dt*(Dzz*Vz2_hat)'*Vz2_hat, ((rvals .* Vr2_hat)')*W1*Vz2_hat);
    [Vr2, S2, Vz2, ~] = truncate_svd(Vr2_hat, S2_hat, Vz2_hat, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vr_dagger3, ~, Vz_dagger3, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, 1.5*dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    Y2 = (((Drr*Vr2*S2*(Vz2')) + (Vr2*S2*((Dzz*Vz2)'))));
    W2 = (Vr0*S0*(Vz0')) + (-1*dt*Y1) + (2*dt*Y2);
      
    % Reduced augmentation
    [Vr_star3, Vz_star3] = reduced_augmentation([Vr_dagger3, Vr2, Vr1, Vr0], [Vz_dagger3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K3 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star3)'*Vz_star3, W2*Vz_star3);
    L3 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*Vr_star3)'*(rvals .* Vr_star3), (W2')*(rvals .* Vr_star3));

    % Get bases
    [Vr_ddagger3, ~] = qr2(K3, rvals); [Vz_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vr3_hat, Vz3_hat] = reduced_augmentation([Vr_ddagger3, Vr2, Vr1, Vr0], [Vz_ddagger3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S3_hat = sylvester(eye(size(Vr3_hat, 2)) - (a_kk*dt*((rvals .* Vr3_hat)')*Drr*Vr3_hat), -a_kk*dt*(Dzz*Vz3_hat)'*Vz3_hat, ((rvals .* Vr3_hat)')*W2*Vz3_hat);
    [Vr3, S3, Vz3, ~] = truncate_svd(Vr3_hat, S3_hat, Vz3_hat, tolerance);

    % Stage 4:
    % Predict V_dagger using B. Euler
    [Vr_dagger4, ~, Vz_dagger4, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, rvals, zvals, Rmat, Zmat, Drr, Dzz, tolerance, rhoM);
    Y3 = (((Drr*Vr3*S3*(Vz3')) + (Vr3*S3*((Dzz*Vz3)'))));
    W3 = (Vr0*S0*(Vz0')) + (-(1/12)*dt*Y1) + ((2/3)*dt*Y2) + (-(1/12)*dt*Y3);
      
    % Reduced augmentation
    [Vr_star4, Vz_star4] = reduced_augmentation([Vr_dagger4, Vr3, Vr2, Vr1, Vr0], [Vz_dagger4, Vz3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K4 = sylvester(eye(size(Drr)) - (a_kk*dt*Drr), -a_kk*dt*(Dzz*Vz_star4)'*Vz_star4, W3*Vz_star4);
    L4 = sylvester(eye(size(Dzz)) - (a_kk*dt*Dzz), -a_kk*dt*(Drr*Vr_star4)'*(rvals .* Vr_star4), (W3')*(rvals .* Vr_star4));

    % Get bases
    [Vr_ddagger4, ~] = qr2(K4, rvals); [Vz_ddagger4, ~] = qr(L4, 0);

    % Reduced Augmentation
    [Vr4_hat, Vz4_hat] = reduced_augmentation([Vr_ddagger4, Vr3, Vr2, Vr1, Vr0], [Vz_ddagger4, Vz3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S4_hat = sylvester(eye(size(Vr4_hat, 2)) - (a_kk*dt*((rvals .* Vr4_hat)')*Drr*Vr4_hat), -a_kk*dt*(Dzz*Vz4_hat)'*Vz4_hat, ((rvals .* Vr4_hat)')*W3*Vz4_hat);
    % [Vr4, S4, Vz4, rank] = truncate_svd(Vr4_hat, S4_hat, Vz4_hat, tolerance);
    [Vr4, S4, Vz4, rank] = LoMaC_mass_only(Vr4_hat, S4_hat, Vz4_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM);
end



function [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance)
    [U, Sigma, V] = svd(S, 0);
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    if (sum(rank) == 0)
        rank = 1;
    end
    Vr = Vr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Vz*V(:, 1:rank);
end

function [Vr, Vz] = reduced_augmentation(Vr_aug, Vz_aug, rvals)
    tolerance = 1e-12;
    [Qr, Rr] = qr2(Vr_aug, rvals);
    [Qz, Rz] = qr(Vz_aug, 0);
    [Ur, Sr, ~] = svd(Rr, 0);
    [Uz, Sz, ~] = svd(Rz, 0);
    rr = find(diag(Sr) > tolerance, 1, 'last');
    rz = find(diag(Sz) > tolerance, 1, 'last');
    rank = max(rr, rz);
    rank = min(rank, min(size(Ur, 2), size(Uz, 2)));
    Vr = Qr*Ur(:, 1:rank);
    Vz = Qz*Uz(:, 1:rank);
end

function [Q, R] = qr2(X, rvals)
    [Q, R] = qr(sqrt(rvals) .* X, 0);
    Q = Q ./ sqrt(rvals);
end

function [U, S, V] = svd2(X, rvals)
    [U, S, V] = svd(sqrt(rvals) .* X, 0);
    U = U./sqrt(rvals);
end

function [Rmat, Zmat, dr, dz] = GetRZ(vmin, vmax, zmin, zmax, Nv, Nz)
    rvals = linspace(vmin, vmax, Nv+1)';
    zvals = linspace(zmin, zmax, Nz+1)';
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    rmid = rvals(1:end-1) + (dr/2);
    zmid = zvals(1:end-1) + (dz/2);
    [Rmat, Zmat] = meshgrid(rmid, zmid);
    Rmat = Rmat';
    Zmat = Zmat';
end



% ------- LoMaC Truncation -------
function [Vr, S, Vz, rank] = LoMaC_mass_only(Vr, S, Vz, R, Z, rvals, zvals, tolerance, rhoM)
    % LoMaC Truncates given maxwellian (assumed Low-Rank) to given tolerance while conserving
    % macroscopic quantities (in this case, only mass for the heat equation).

    Nr = numel(rvals); Nz = numel(zvals);
    dvr = rvals(2) - rvals(1);
    dvz = zvals(2) - zvals(1);

    % Step 1: Integrate to calculate macro quantities
    p = 2*pi*dvr*dvz*sum(sum(((Vr) * S * (Vz)') .* R));

    % Step 2: Scale by maxwellian to ensure inner product is well defined
    % (f -> 0 as v -> infinity)
    dropoff = 10;
    wr = exp(-dropoff*(rvals.^2));
    wz = exp(-dropoff*(zvals.^2));    

    % Step 3: Orthogonal projection
    % bases: 1, v, v.^2 - c
    % c = (dvr*sum(rvals.^2.*wr.*rvals))/(dvr*sum(wr.*rvals)) + (dvz*sum(zvals.^2.*wz))/(dvz*sum(wz));

    w_norm_1_squared = 2*pi*dvr*dvz*sum(rvals .* wr)*sum(wz);
    % w_norm_v2_squared = (2*pi*dvr*dvz*sum(sum((R.^2 + Z.^2 - c).^2 .* exp(-R.^2 - Z.^2) .* R)));
    
    f1_proj_S_mtx11 = (p / w_norm_1_squared);

    proj_basis_r = wr.*[ones(Nr, 1)];
    proj_basis_z = wz.*[ones(Nz, 1)];
    f1_proj_S_mtx   = f1_proj_S_mtx11;

    % f2 = f - f1 (do it via SVD)
    f2_U = [Vr, proj_basis_r];
    f2_S = blkdiag(S, -f1_proj_S_mtx);
    f2_V = [Vz, proj_basis_z];

    % QR factorize
    [f2_Vr, f2_S, f2_Vz, ~] = truncate(f2_U, f2_S, f2_V, rvals, tolerance);

    f2 = f2_Vr*f2_S*f2_Vz';

    % compute Pn(Te(f)) to ensure moments are kept
    trun_f2_proj_S_mtx11 = 2*pi*dvr*dvz*(sum(sum(f2.*R)) / w_norm_1_squared);% - (c*sum(sum(f2.*(R.^2 + Z.^2 - c) .* R)) / w_norm_v2_squared));
    trun_f2_proj_S_mtx = trun_f2_proj_S_mtx11;

    % compute fM
    fM_proj_S_mtx11 = (rhoM ./ w_norm_1_squared);
    fM_proj_S_mtx = fM_proj_S_mtx11;

    f_mass_S = fM_proj_S_mtx - trun_f2_proj_S_mtx;

    [Vr, S, Vz, rank] = truncate([proj_basis_r, f2_Vr], blkdiag(f_mass_S, f2_S), [proj_basis_z, f2_Vz], rvals, 1e-14);
end

function [Vr, S, Vz, rank] = truncate(Vr_aug, S_aug, Vz_aug, rvals, tolerance)
    [Qr, Rr] = qr2(Vr_aug, rvals); [Qz, Rz] = qr(Vz_aug, 0);
    [U, Sigma, V] = svd(Rr*S_aug*Rz', 0); 
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    if numel(rank) == 0
        rank = 1;
    end
    Vr = Qr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Qz*V(:, 1:rank);
end









% %%
% % l1 decay
% 
% semilogy(tvals, l1, 'LineWidth', 1.5); hold on;
% % xlabel('t'); ylabel('L_1(f(V_r, V_z))'); title('L_1 drive to equilibrium solution');
% 
% % Relative entropy
% semilogy(tvals, relative_entropy,  'LineWidth', 1.5);
% % xlabel('t'); ylabel('Relative entropy'); title('Relative entropy decay');
% hold off
% legend('L1 decay', 'Relative entropy')
% xlabel('Time')
% ylabel('Magnitude')
% ylim([1e-15, 1e0]);
% title('L_1 drive to equilibrium solution and relative entropy decay');
% return
% %%
% % Positivity
% figure(4); clf; plot(tvals, min_vals, 'green-', 'LineWidth', 1.5);
% xlabel('t'); ylabel('min(f(V_r, V_z))'); title('Minimum values of numerical solution over time');
% 
% % Mass
% figure(6); clf; plot(tvals(2:end), abs(mass(2:end)-mass(1))/mass(1), 'red-', 'LineWidth', 1.5);
% xlabel('t'); ylabel('relative mass'); title('Relative mass of numerical solution over time');
% figure(7); clf; plot(tvals(2:end), abs(Jzvals(2:end)-Jzvals(1)), 'red-', 'LineWidth', 1.5);
% xlabel('t'); ylabel('Absolute error (Uz)'); title('Absolute error of bulk velocity over time');
% figure(8); clf; plot(tvals(2:end), abs(E(2:end)-E(1))/E(1), 'red-', 'LineWidth', 1.5);
% xlabel('t'); ylabel('Relative error (Energy)'); title('Relative energy of numerical solution over time');
% 
% 
% 
% % Rank plot
% figure(9); clf;
% plot(tvals, ranks, 'black-', 'LineWidth', 1.5); hold on;
% % plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
% % plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
% xlabel('time'); ylabel('rank'); title('Rank plot over time');
% legend('Backward Euler', 'RK2', 'RK3');
