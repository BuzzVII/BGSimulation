set(0,'DefaultAxesFontSize',20) 
subplot(2,1,1)
t = 1/10000*(1:length(rates_full));
plot(t,abs(rates_0 - rates_full)./rates_full,'-.','linewidth', 2)
hold on; plot(t,abs(rates_1 - rates_full)./rates_full,'g','linewidth', 2)
xlabel('time (s)')
ylabel('fraction error')
axis([0 0.5 0 2e-3])
subplot(2,1,2)
t = 1/10000*(1:length(rates_full_lin));
plot(t,abs(rates_0_lin - rates_full_lin)./rates_full_lin,'-.','linewidth', 2)
hold on; plot(t,abs(rates_1_lin - rates_full_lin)./rates_full_lin,'g','linewidth', 2)
xlabel('time (s)')
ylabel('fraction error')
axis([0 5 0 1e-4])
legend('0th order','1st order')