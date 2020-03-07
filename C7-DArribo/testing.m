t = 1:100;
x = e.^(-t/10);
y = [zeros(1,15) x(1:end-15)];
plot(t,x)