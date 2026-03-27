clear; clc;

%% PARAMETERS
N = 5000;
t_values = linspace(-sqrt(5)/(2*sqrt(2)), sqrt(5)/(2*sqrt(2)), N);

states = zeros(4,4,N);
labels = zeros(N,1);   % optional PPT label

for k = 1:N
    
    t = t_values(k);
    
    % ===== DEFINE YOUR MATRIX HERE =====
    % Replace this with your actual formula
    
    rho = (1/2) *[ 5/4,  0,     0,     t;
                0,     0,  0,     0;
                0,     0,     1/4,  0;
                t,     0,     0,     1/2 ];
    % ==================================
    
    states(:,:,k) = rho;
    
    % Optional: PPT label
    pt = PartialTranspose(rho);
    if min(eig(pt)) < -1e-10
        labels(k) = 1;   % entangled
    else
        labels(k) = 0;   % separable
    end
end
save('rudolph_twoqubit_dataset.mat', ...
     'states','labels','t_values','-v7.3');

disp('Dataset created successfully.');
