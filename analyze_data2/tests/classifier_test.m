


load temp

x=y;

dx=diff(x);
dy=diff(y);


figure; plot([x]); xlabel('time');
hold on; plot(dx,'k'); legend('x','dx/dt')
 
figure; hist(x,100); xlabel('x bin'); ylabel('N'); title('Histogram of dx/dt')
figure; hist(dx,100);  xlabel('x bin'); ylabel('N'); title('Histogram of dx/dt')



