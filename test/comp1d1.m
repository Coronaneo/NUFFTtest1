Nlist = 2.^(7:12);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timNUFFTAppnyu=zeros(size(Nlist));
timeM=zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));
tol=1e-12;
num = 2;

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x=sort(N*rand(N,1));
    x1=2*pi*x/N;
    
    tic;
    for cnt = 1:num
    nufftfun = nufftI(x,-1,N,15,tol);
    end
    timNUFFTFact(it) = toc/num;
    
    c=rand(N,1);
    
    tic;
    for cnt = 1:num
    nufftc = nufftfun(c);
    end
    timNUFFTApp(it) = toc/num;
    
    tic;
    for cnt = 1:num
    %fk=nufft1d1(N,x1,c,-1,tol,N)*N;
    %fk=fftshift(fk);
    end
    timNUFFTAppnyu(it)=toc/num;
    
    
    k=-N/2:(N/2-1);
    k = k(:);
    [fhatM,ffun] = DeCom_NUFFT1D_I(c,x,k,tol);
    tic;
    for cnt = 1:num
        fhatM = ffun(c);
    end
    timeM(it) = toc/num;
    
    
    %errNUFFT(it)=norm(nufftc(1:ceil(N/2))-fk(1:ceil(N/2)),2)/norm(nufftc(1:ceil(N/2)),2);
    
end

timNUFFTApp
timNUFFTAppnyu
timeM
timecomp=timeM./timNUFFTAppnyu
fid=fopen('./result1d1/time1d1YH.mat','at');
fprintf(fid,'% -f\n',timeM);
fclose(fid);
fid=fopen('./result1d1/time1d1NYU.mat','at');
fprintf(fid,'% -f\n',timNUFFTAppnyu);
fclose(fid);
fid=fopen('./result1d1/time1d1LR.mat','at');
fprintf(fid,'% -f\n',timNUFFTApp);
fclose(fid);
