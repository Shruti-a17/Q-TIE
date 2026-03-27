clear; clc;

%% ==============================
% PARAMETERS
%% ==============================

r = 3;                 % rank (1,2,3,4)
N_ent = 6000;          % total entangled states

nbins = 60;            
edges = linspace(0,1,nbins+1);
per_bin = N_ent/nbins;

labels = ones(N_ent,1);   % all entangled

%% ==============================
% PARALLEL SETUP
%% ==============================

delete(gcp('nocreate'));   % reset pool (important)
parpool;

%% ==============================
% STORAGE (CELL FOR PARFOR)
%% ==============================

states_cell = cell(nbins,1);
conc_cell   = cell(nbins,1);

%% ==============================
% MAIN PARALLEL LOOP
%% ==============================

parfor b = 1:nbins

    Cmin = edges(b);
    Cmax = edges(b+1);

    fprintf('Bin %d: (%.3f, %.3f)\n', b, Cmin, Cmax);

    local_states = zeros(4,4,per_bin);
    local_conc   = zeros(per_bin,1);

    count = 0;

    while count < per_bin

        % ---- Generate random rank-r state ----
        rho = random_rank_r_state(r);

        % ---- Rank check ----
        if rank(rho) ~= r
            continue
        end

        % ---- Concurrence ----
        C = Concurrence(rho);

        if C <= Cmin || C >= Cmax
            continue
        end

        % ---- PPT check (entangled) ----
        pt = PartialTranspose(rho);

        if min(eig(pt)) < -1e-10

            count = count + 1;
            local_states(:,:,count) = rho;
            local_conc(count) = C;

        end
    end

    states_cell{b} = local_states;
    conc_cell{b}   = local_conc;

end

%% ==============================
% MERGE DATA
%% ==============================

states = cat(3, states_cell{:});
conc   = vertcat(conc_cell{:});

%% ==============================
% VERIFY
%% ==============================

disp('Final dataset size:')
disp(size(states))   % should be 4 x 4 x 6000
disp(length(conc))

%% ==============================
% SAVE DATASET
%% ==============================

filename = sprintf('rank_%3_entangled_6000.mat', r);

save(filename, 'states', 'conc', 'labels', '-v7.3');

disp(['Dataset saved as ', filename]);

%% ==============================
% OPTIONAL: CHECK DISTRIBUTION
%% ==============================

figure;
histogram(conc,50);
xlabel('Concurrence');
ylabel('Count');
title('Concurrence Distribution');
grid on;

%% ==============================
% FUNCTION: RANDOM RANK-r STATE
%% ==============================

function rho = random_rank_r_state(r)

    % ---- Bias eigenvalues (critical for speed) ----
    x = rand(r,1).^3;   % helps high concurrence bins
    x = x / sum(x);

    D = zeros(4);
    D(1:r,1:r) = diag(x);

    % ---- Haar random unitary ----
    X = randn(4) + 1i*randn(4);
    [Q,R] = qr(X);
    d = diag(R);
    U = Q * diag(d./abs(d));

    rho = U * D * U';

end
