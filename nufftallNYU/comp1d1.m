Nlist = 2.^(7:20);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timNUFFTAppnyu=zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));
eps=1e-12;

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x=sort(N*rand(N,1));
    x1=2*pi*x/N;
    
    tic;
    nufftfun = nufftI(x,-1,N);
    timNUFFTFact(it) = toc;
    
    c=rand(N,1);
    
    tic;
    nufftc = nufftfun(c);
    timNUFFTApp(it) = toc;
    
    tic;
    fk=nufft1d1(N,x1,c,-1,eps,N)*N;
    fk=fftshift(fk);
    timNUFFTAppnyu(it)=toc;
    
    errNUFFT(it)=norm(nufftc(1:ceil(N/2))-fk(1:ceil(N/2)),2)/norm(nufftc(1:ceil(N/2)),2);
    
end

sub=timNUFFTAppnyu./timNUFFTApp;
figure
loglog(Nlist,sub,'b');
hold on
loglog(Nlist,ones(size(Nlist)),'r');
xlabel('N'),ylabel('NYUtime/Ourtime'),title('One-dimensional-TypeI Comparison ')
errNUFFT