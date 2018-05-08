clear;

%MER paramters
N = 1;
T = 10.0;
fs = 24000.0;
dt = 1/1000;
da = 1.0;

omegas = 10.^(-2:0.1:2);
rate_all = {};
shape_all = {};
cov_all = {};

%Neuron parameters
It = dlmread('apcurrent24k.dat');
It = It / max(abs(It));

sim_N = 300;
for omega = omegas
    cov = zeros(1,sim_N);
    shape = zeros(1,sim_N);
    rate = zeros(1,sim_N);
    %BG parameters
    v0 = 14.0;
    v_e = @(t) omega*2*pi*cos(2*pi*20*t); %cortical rate of change
    W = @(t) 0;%omega*randn(); %cortical noise function
    c = 0; %1/correlation time of noise
    %U =  [9.293; -18.4203; 3.857; -4.0643; 16.2714; -194.239; 11.39; -135.967; -1.544; -18.446; 14.0; 14.0; 14.0];
    U = [9.293; 0; 3.857; 0; 14.0; 0; 9.8; 0; -1.78; 0; v0; v0; v0];
    for ind = 1:sim_N
        if mod(ind,10) == 0
            fprintf('%d out of %d for omega %.2f',ind,sim_N, omega);
        end
        [rates, sol, t] = bg_sim(v_e, W, c, U, T, dt, da);

        %the simulation
        stn_rate = interp1(t, rates(:,5), (0:1/fs:T)); 
        neuron_superposition = generate_timing(T, fs, stn_rate, N);
        mer = conv(neuron_superposition, It);
        %mer = mer + 0.01*randn(size(mer));
        mer = mer - mean(mer);
        h = fit_weibull(mer, fs, It, 5);
        shape(ind) = h(2);
        rate(ind) = h(1);
        isi = diff(find(neuron_superposition >= 1));
        cov(ind) =  std(isi)/mean(isi);
    end
    shape_all{end+1} = shape;
    rate_all{end+1} = rate;
    cov_all{end+1} = cov;
end
labels = {'Correlation Parameter'};
plot_simulations(shape_all, omegas, labels)