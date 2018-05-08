load('weiner_sim.mat')
labels = {'Diffusion parameter'};
plot_simulations(shape_all(end-18:end), omegas(end-18:end), labels)
load('correlated_sim.mat')
labels = {'Correlation parameter'};
plot_simulations(shape_all, omegas, labels)
% load('wsin(t).mat')
% labels = {'Amplitude parameter'};
% plot_simulations(shape_all, omegas, labels)
load('sin(wt).mat')
labels = {'Frequency parameter'};
plot_simulations(shape_all, omegas, labels)