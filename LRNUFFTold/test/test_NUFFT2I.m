Nlist = 2.^(4:8);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timFFTApp = zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x = N*rand(N^2,2);
    [k1,k2] = ndgrid(0:N-1);
    k = [k1(:) k2(:)];
    
    tic;
    nufft2fun = nufft2I(x,-1,N);
    timNUFFTFact(it) = toc;
    
    idx = randsample(N^2,128);
    
    c = randn(N^2,2);
    
    tic;
    nufft2c = nufft2fun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fftshift(fft2(reshape(c,N,N,[])));
    timFFTApp(it) = toc;
    
    
    nuFFTDensefun = @(k,x)exp(-1i*2*pi*k*x'/N);
    nuFFTMat = nuFFTDensefun(k(idx,:),x);
    errNUFFT(it) = norm(nufft2c(idx,:)-nuFFTMat*c)/norm(nuFFTMat*c);
    
end

figure
loglog(Nlist,timNUFFTFact);

figure
loglog(Nlist,timNUFFTApp);
hold on;
loglog(Nlist,timFFTApp);
