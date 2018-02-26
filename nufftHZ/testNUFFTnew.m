function testNUFFTnew()
% This code tests the idea of the NUFFT
close all;
clear all;

N = 2^10;
tol = 1e-13;
eps = 0.01;%1/N/4;
x = (0:1/N:(1-1/N))+rand(1,N)*eps;
x = mod(x,1);
x = sort(x);
x = x(:);
%x = sort(rand(1,N));
k=-N/2:(N/2-1);
%k = 0:N-1;
k = k(:);
A = exp(-2*pi*i*k*x');
F = 100;
f = cos(2*pi*F*x);
ff = cos(2*pi*F*(0:1/N:(1-1/N)));
ffhat = fft(ff(:));
num = 10;

fhat = ffhat;
fhat = A*f;
% Mine
[fhatM,ffun] = DeCom_NUFFT1D_I(f,x,k,tol*10);
tic;
for cnt = 1:num
    fhatM = ffun(f);
end
timeM = toc;
timeM = timeM/num;
% Alex's
tic;
[fhatA,p] = nufft(f,x*N,1);
timeA = toc;

% Greengard's
tic;
for cnt = 1:num
    fhatG = nufft1d1(N,x*2*pi,f,-1,tol,N);
    fhatG = fhatG*(N);
   % fhatG = fftshift(fhatG);
end
timeG = toc;
timeG = timeG/num;

errM = norm(fhatM-fhat)/norm(fhat);
errA = norm(fhatA-fhat)/norm(fhat);
errG = norm(fhatG-fhat)/norm(fhat);
fprintf('error of Greengard is %.3e\n',errG);
fprintf('error of Alex is %.3e\n',errA);
fprintf('error of mine is %.3e\n\n',errM);

errM = norm(fhatM-ffhat)/norm(ffhat);
errA = norm(fhatA-ffhat)/norm(ffhat);
errG = norm(fhatG-ffhat)/norm(ffhat);
fprintf('error of Greengard is %.3e\n',errG);
fprintf('error of Alex is %.3e\n',errA);
fprintf('error of mine is %.3e\n\n',errM);

fprintf('time of Greengard is %.3e\n',timeG);
fprintf('time of Alex is %.3e\n',timeA);
fprintf('time of mine is %.3e\n',timeM);


figure;subplot(5,2,1);plot(real(fhat));subplot(5,2,2);plot(imag(fhat));
subplot(5,2,3);plot(real(fhatA));subplot(5,2,4);plot(imag(fhatA));title('Alex');
subplot(5,2,5);plot(real(fhatG));subplot(5,2,6);plot(imag(fhatG));title('Greengard');
subplot(5,2,7);plot(real(fhatM));subplot(5,2,8);plot(imag(fhatM));title('Mine');
subplot(5,2,9);plot(real(ffhat));subplot(5,2,10);plot(imag(ffhat));title('Uniform');
timeG/timeM
end

function f = nudft1( c, omega)

f = zeros(size(omega,1),1);
for j = 1:numel(f)
    f(j) = exp(-2*pi*1i*(j-1)/size(omega,1)*omega.')*c;
end
end
