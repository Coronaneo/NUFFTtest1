Nlist = 2.^(7:15);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timFFTApp = zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x = (N-1)*rand(N,1);
    k = (N-1)*rand(N,1);
    
    tic;
    nufftfun = nufftIII(k,x,-1,N);
    timNUFFTFact(it) = toc;
    
    idx = randsample(N,128);
    
    c = randn(N,1);
    
    tic;
    nufftc = nufftfun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fftshift(fft(c));
    timFFTApp(it) = toc;
    
    nuFFTDensefun = @(k,x)exp(-1i*2*pi*k*x'/N);
    nuFFTMat = nuFFTDensefun(k(idx,:),x);
    errNUFFT(it) = norm(nufftc(idx)-nuFFTMat*c)/norm(nuFFTMat*c);
    
end

figure
loglog(Nlist,timNUFFTFact);

figure
loglog(Nlist,timNUFFTApp);
hold on;
loglog(Nlist,timFFTApp);
