Nlist = 2.^(7:12);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timNUFFTAppold = zeros(size(Nlist));
timNUFFTAppnyu=zeros(size(Nlist));
timeM=zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));
errNUFFT1 = zeros(size(Nlist));
errNUFFT2 = zeros(size(Nlist));
tol=1e-6;
num = 100;
iflag=-1;

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x=sort(N*rand(N,1));
    x1=2*pi*x/N;
    
    tic;
    
    nufftfun = nufftII(x,iflag,N,13,tol);
    
    timNUFFTFact(it) = toc;
    
    c=rand(N,1);
    
    tic;
    for cnt = 1:num
    nufftc = nufftfun(c);
    end
    timNUFFTApp(it) = toc/num;
    
    %for cnt = 1:num
    %nufftfun = nufftIIold(x,iflag,N,15,tol);
    %end
    %timNUFFTFact(it) = toc/num;
    
    %tic;
    %for cnt = 1:num
    %nufftc1 = nufftfun(c);
    %end
    %timNUFFTAppold(it) = toc/num;
    
    tic;
    for cnt = 1:num
    cj1 = nufft1d2(N,x1,iflag,tol,N,c);
    
    end
    timNUFFTAppnyu(it)=toc/num;
    
    
    k=-N/2:(N/2-1);
    [fhatM,ffun] = DeCom_NUFFT1D_II(c,x/N,k,tol);
    
    tic;
    %profile on
    for cnt = 1:num
        fhatM = ffun(c);
    end
    %profile report
    timeM(it) = toc/num;
    
    
    errNUFFT(it)=norm(cj1 - fhatM,2)/norm(fhatM,2);
    %errNUFFT1(it)=norm(nufftc-cj1,2)/norm(nufftc,2);
    errNUFFT1(it)=norm(nufftc - fhatM,2)/norm(nufftc,2);
end

timNUFFTApp
%timNUFFTAppold
%timNUFFTAppnyu
timeM
errNUFFT,
errNUFFT1
figure
loglog(Nlist,timNUFFTApp);
hold on;
loglog(Nlist,timeM);
xlabel('N'),ylabel('time : s'),title('Comparison,1D type II,tol =1e-6'),legend('timeYZ','timeHZ','Location','northwest')
%errNUFFT1
%errNUFFT2
%timecomp=timeM./timNUFFTAppnyu
fid=fopen('./result1d1/time1d1YH.mat','at');
fprintf(fid,'% -f\n',timeM);
fclose(fid);
%fid=fopen('./result1d1/time1d1LRold.mat','at');
%fprintf(fid,'% -f\n',timNUFFTAppold);
%fclose(fid);
%fid=fopen('./result1d1/time1d1NYU.mat','at');
%fprintf(fid,'% -f\n',timNUFFTAppnyu);
%fclose(fid);
fid=fopen('./result1d1/time1d1LR.mat','at');
fprintf(fid,'% -f\n',timNUFFTApp);
fclose(fid);