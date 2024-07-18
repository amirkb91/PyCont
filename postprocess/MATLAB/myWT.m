% function [finst,wtinst,time,y,freq] = myWT(x,fs,fi,ff,nf,f0,pad)
%
% Morlet-based wavelet transform computation 
%
% Inputs:
%    x        : signal to be processed
%    fs       : sampling frequency
%    fi/ff    : lower/upper bounds of the working frequency interval
%    nf       : number of frequency steps
%    f0       : mother wavelet frequency
%    pad      : zero-padding parameter (default = 0).
%
% Outputs:
%    finst    : instantaneous frequency vector
%    wtinst   : wavelet transform evaluated at finst (i.e. envelope of the wt)
%    time     : time vector
%    y        : wavelet transform
%    freq     : output frequency vector
%
% 
%
function [finst,wtinst,time,y,freq] = myWT(x,fs,fi,ff,nf,f0,pad)

if nargin==7, pad=0; end

h = 1/fs;
 
hf = (ff-fi)/nf;
freq = fi:hf:ff;

a = zeros(nf+1,1);
for i = 0:nf
  a(i+1) =  f0/(fi+i*hf);
end

Na = length(a)-1;

k = 2^pad;
NX = length(x);
NX2 = nextpow2(NX);
N = 2^NX2;
N = k*N;

time = 0:h:(N-1)*h;

f = 0:fs/(N-2):fs/2;
omega = 2*pi*f;

filter = zeros(Na+1,N/2);
for j = 1:(Na+1)   
   int = a(j).*omega;
   for jj = 1:N/2
      if int(jj) == 0
            filter(j,jj) = 0;
      else
            filter(j,jj) = sqrt(2*a(j))*exp(-0.5*(int(jj)-2*pi*f0).^2);
      end
   end
end

X = fft(x.',N);

y = zeros(N,Na+1);
for j = 1:(Na+1)
     ker(1:N/2) = conj(filter(j,1:N/2)).*X(1:N/2);
     iker = ifft(ker,N);
     y(1:N,j) = (iker(1:N)).';
     
end 
 
mod = abs(y);
% ph = unwrap(angle(wt));

[wtinst imax] = max(mod,[],2);

finst = f0./a(imax);

return