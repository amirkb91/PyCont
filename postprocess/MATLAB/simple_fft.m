close all; clear all;

% Datafile to read from
file = "Node16_Y_lie23.dat";
data = load(file);

t = data(:,1);
x = data(:,2);
figure; plot(t,x,'.-b','LineWidth', 1.2); grid

% FFT
figure;
Fs = 1/(t(2)-t(1));
L = length(t)-1;
L2 = floor(L/2);
ZF = fft(x);
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;

% first point is steady state gain
plot(f(2:end),P(2:end),'b','LineWidth', 1.2); grid;
[pk_P, pk_ind] = findpeaks(P);
pk_f = f(pk_ind);
xlabel('f (Hz)'); ylabel('|Y(f)|'); title('Single-Sided Amplitude Spectrum of Y(t)');
