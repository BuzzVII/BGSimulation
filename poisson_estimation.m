N = 100;
range = 100;
scale = 1;
t = 10;

means = zeros(range,1);
stds = zeros(range,1);
times = zeros(range,1);
counts = zeros(range,1);
mean_stds = zeros(range,1);
std_stds = zeros(range,1);
time_stds = zeros(range,1);
counts_stds = zeros(range,1);
for ind = 1:range
    exp_mean = zeros(N,1);
    exp_std = zeros(N,1);
    exp_time = zeros(N,1);
    exp_counts = zeros(N,1);
    for j = 1:N
        x = randraw('weibull', [0,0.5,1/gamma(1+1/0.5)], scale*ind);
        X = cumsum(x);
        exp_counts(j) = sum(X<t);
        exp_mean(j) = mean(x);
        exp_std(j) = std(x);
        exp_time(j) = sum(x);
    end
    means(ind) = mean(exp_mean);
    stds(ind) = mean(exp_std);
    times(ind) = mean(exp_time);
    counts(ind) = mean(exp_counts);
    mean_stds(ind) = std(exp_mean);
    std_stds(ind) = std(exp_std);
    time_stds(ind) = std(exp_time);
    counts_stds(ind) = std(exp_counts);
end

close all
set(0,'DefaultAxesFontSize',20)
set(0,'defaultlinelinewidth',2)
figure()
loglog(means)
xlabel('samples')
ylabel('\mu')
figure()
loglog(stds)
xlabel('samples')
ylabel('\sigma')
figure()
loglog(times)
xlabel('samples')
ylabel('time \mu')
figure()
loglog(mean_stds)
xlabel('samples')
ylabel('\mu \sigma')
figure()
loglog(std_stds)
xlabel('samples')
ylabel('\sigma \sigma')
figure()
loglog(time_stds)
xlabel('samples')
ylabel('time \sigma')

figure()
plot(mean_stds./means)
xlabel('samples')
ylabel('\mu % accuracy')

figure()
loglog(counts)
xlabel('samples')
ylabel('counts')
figure()
loglog(mean_stds)
xlabel('samples')
ylabel('counts \sigma')