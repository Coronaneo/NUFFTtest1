% This code tests the idea of the sparse and partial NUFFT
clear all;
close all;

N = 2^15;  %18    %% sampling size
opt.bandWidth = 200;
opt.Ls = 4000;
k = -N/2:(N/2-1);
tol = 1e-13;
x0 = 0:1/opt.Ls:(1-1/opt.Ls);

x = [0:N-1]/N;
L = 150;
sh = @(x) gen_shape2(x,3);
figure;subplot(2,2,1);plot(x0,sh(x0));
amp = 0.01;
xx = (x + amp*sin(2*pi*x));
xx = xx - xx(1);
fq = (1+amp*2*pi*cos(2*pi*x))*L;
sig = sh(L*xx);

%**********************************************************
tic;
fphah=nufft1d1(N,xx*2*pi,sig,-1,1e-14,N);
% [fphah,ffun] = DeCom_NUFFT1D_I(sig(:),xx(:),k(:),tol*10);
% fphah = fphah.'/N;
fphah = transpose(fphah)*sqrt(N);
vector1 = (0:L:N/2-L/2)';
vector1 = vector1(1:min(opt.bandWidth+1,numel(vector1)));
vector2 = (-L:-L:-N/2+L/2)';
vector2 = vector2(end:-1:1);
vector2 = vector2(1:min(opt.bandWidth,numel(vector2)));
vector=[vector2;vector1];
K=numel(vector);
u = fphah(vector+floor(N/2)+1);
shapehat =zeros(opt.Ls,1);
if opt.Ls>=K
    st = opt.Ls/2+1-floor(K/2);
    shapehat(st:st+K-1)=u;
else
    st = floor(K/2)+1-floor(opt.Ls/2);
    shapehat=u(st:st+opt.Ls-1);
end
% scaling s.t. shape is independent of opt.Ls
shapeTemp = ifft(fftshift(shapehat))*opt.Ls/sqrt(N);
shape=real(shapeTemp)';
shift = mean(shape);
shape = shape - shift;
timeG = toc;

subplot(2,2,2);plot(x0,shape);

%**********************************************************
tic;
kk = 0:L:(N/2-1); tmp = -L:-L:-N/2; kk = [tmp,kk];
Nk = numel(kk);
st = randi([1,N-Nk+1],1,1);
st = floor(st/Nk)*Nk;

fphah=nufft1d1(Nk,xx(st:st+Nk-1)*2*pi,sig(st:st+Nk-1),-1,1e-14,Nk);
u = transpose(fphah)*sqrt(Nk);
shapehat =zeros(opt.Ls,1);
if opt.Ls>=Nk
    st = opt.Ls/2+1-floor(Nk/2);
    shapehat(st:st+Nk-1)=u;
else
    st = floor(Nk/2)+1-floor(opt.Ls/2);
    shapehat=u(st:st+opt.Ls-1);
end
% scaling s.t. shape is independent of opt.Ls
shapeTemp = ifft(fftshift(shapehat))*opt.Ls/sqrt(N);
shape=real(shapeTemp)';
shift = mean(shape);
shape = shape - shift;
timeG = toc;

subplot(2,2,3);plot(x0,shape);

%**********************************************************

[fhat,ffun] = DeCom_SPNUFFT1D_I(sig(:),xx(:),k(:),opt.bandWidth,L,tol);
tic;
u = ffun(sig(:));
K = numel(u);
shapehat =zeros(opt.Ls,1);
if opt.Ls>=K
    st = opt.Ls/2+1-floor(K/2);
    shapehat(st:st+K-1)=u/K;
else
    st = floor(K/2)+1-floor(opt.Ls/2);
    shapehat=u(st:st+opt.Ls-1)/K;
end
shapeTemp = ifft(fftshift(shapehat))*opt.Ls/sqrt(N);
shape=real(shapeTemp)';
shift = mean(shape);
shape = shape - shift;
timeM = toc;

subplot(2,2,4);plot(x0,shape);
timeG/timeM
