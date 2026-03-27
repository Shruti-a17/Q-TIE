clear; clc;

%% ===============================
% PARAMETERS
%% ===============================

Total = 25000;
Nbins = 10;
N_per_bin = Total / Nbins;

P_min = 0.25;
P_max = 1.0;

purity_edges = linspace(P_min, P_max, Nbins+1);

states = zeros(4,4,Total);
labels = zeros(Total,1);  % all PPT ? label 0
purities = zeros(Total,1);

index = 1;

%% ===============================
% GENERATE STATES PER PURITY BIN
%% ===============================

for b = 1:Nbins
    
    P_low  = purity_edges(b);
    P_high = purity_edges(b+1);
    
    fprintf('Generating purity bin (%.3f, %.3f)\n', P_low, P_high);
    
    count = 0;
    
    while count < N_per_bin
        
        % --- Random separable state ---
        % Mixture of 4 random product states
        
        rho = zeros(4);
        probs = rand(4,1);
        probs = probs / sum(probs);
        
        for j = 1:4
            psiA = randn(2,1) + 1i*randn(2,1);
            psiA = psiA / norm(psiA);
            
            psiB = randn(2,1) + 1i*randn(2,1);
            psiB = psiB / norm(psiB);
            
            psi = kron(psiA, psiB);
            rho = rho + probs(j)*(psi*psi');
        end
        
        P = real(trace(rho*rho));
        
        if P >= P_low && P < P_high
            
            states(:,:,index) = rho;
            purities(index) = P;
            labels(index) = 0;  % PPT state
            
            index = index + 1;
            count = count + 1;
        end
    end
    
    fprintf('Finished bin %d\n', b);
end

disp('Generation complete.');

%% ===============================
% SHUFFLE
%% ===============================

perm = randperm(Total);
states = states(:,:,perm);
labels = labels(perm);
purities = purities(perm);

%% ===============================
% SAVE
%% ===============================

dataset.states = states;
dataset.labels = labels;
dataset.purity = purities; 
dataset.N = Total;

save('two_qubit_25k_PPT_uniform_purity.mat','dataset','-v7.3');

disp('Dataset saved successfully.');
