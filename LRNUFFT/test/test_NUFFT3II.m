Nlist = 2.^(3:5);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timFFTApp = zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x = N*rand(N^3,3);
    [k1,k2,k3] = ndgrid(0:N-1,0:N-1,0:N-1);
    k = [k1(:) k2(:) k3(:)];
    
    tic;
    nufft3fun = nufft3II(x,1,N);
    timNUFFTFact(it) = toc;
    
    idx = randsample(N^3,128);
    
    c = randn(N^3,5);
    tic;
    nufft3c = nufft3fun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fftshift(ifft3(reshape(c,N,N,N,[])));
    timFFTApp(it) = toc;
    
    
    nuFFTDensefun = @(x,k)1/N^3*exp(1i*2*pi*x*k'/N);
    nuFFTMat = nuFFTDensefun(x(idx,:),k);
    errNUFFT(it) = norm(nufft3c(idx,:)-nuFFTMat*c)/norm(nuFFTMat*c);
    
end

figure
loglog(Nlist,timNUFFTFact);

figure
loglog(Nlist,errNUFFT);

figure
loglog(Nlist,timNUFFTApp,'b');
hold on;
loglog(Nlist,timFFTApp,'r');