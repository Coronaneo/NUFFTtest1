Nlist = 2.^(3:4);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timFFTApp = zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x = (N-1)*rand(N^3,3);
    k = (N-1)*rand(N^3,3);
    
    tic;
    nufft3fun = nufft3III(k,x,-1,N);
    timNUFFTFact(it) = toc;
    
    idx = randsample(N^3,128);
    
    c = randn(N^3,2);
    
    tic;
    nufft3c = nufft3fun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fftshift(fft3(reshape(c,N,N,N,[])));
    timFFTApp(it) = toc;
    
    nuFFTDensefun = @(k,x)exp(-1i*2*pi*k*x'/N);
    nuFFTMat = nuFFTDensefun(k(idx,:),x);
    errNUFFT(it) = norm(nufft3c(idx,:)-nuFFTMat*c)/norm(nuFFTMat*c);
    
end

figure
loglog(Nlist,timNUFFTFact);

figure
loglog(Nlist,timNUFFTApp);
hold on;
loglog(Nlist,timFFTApp);
