Nlist = 2.^(7:15);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timFFTApp = zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x = N*rand(N,1);
    k = (0:N-1)';
    
    tic;
    nufftfun = nufftII(x,1,N);
    timNUFFTFact(it) = toc;
    
    idx = randsample(N,128);
    
    c = randn(N,5);
    
    tic;
    nufftc = nufftfun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fftshift(ifft(c));
    timFFTApp(it) = toc;
    
    
    nuFFTDensefun = @(x,k)1/N*exp(1i*2*pi*x*k'/N);
    nuFFTMat = nuFFTDensefun(x(idx,:),k);
    errNUFFT(it) = norm(nufftc(idx,:)-nuFFTMat*c)/norm(nuFFTMat*c);
    
end

figure
loglog(Nlist,timNUFFTFact);

figure
loglog(Nlist,timNUFFTApp);
hold on;
loglog(Nlist,timFFTApp);
