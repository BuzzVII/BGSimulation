set(0,'DefaultAxesFontSize',20) 
subplot(2,2,1)
plot(t,sinwt(:,5),'linewidth',2)
xlabel('time (s)')
ylabel('rate (s^{-1})')
axis([1 1.5 20 60])
subplot(2,2,2)
plot(t,wsint(:,5),'linewidth',2)
xlabel('time (s)')
ylabel('rate (s^{-1})')
axis([1 1.5 20 70])
subplot(2,2,3)
plot(t,weiner(:,5),'linewidth',2)
xlabel('time (s)')
ylabel('rate (s^{-1})')
axis([0 5 20 90])
subplot(2,2,4)
plot(t,correlated(:,5),'linewidth',2)
xlabel('time (s)')
ylabel('rate (s^{-1})')
axis([0 5 20 70])