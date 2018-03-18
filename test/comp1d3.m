Nlist = 2.^(7:16);

timNUFFTFact = zeros(size(Nlist));
timNUFFTApp = zeros(size(Nlist));
timNUFFTAppold = zeros(size(Nlist));
timNUFFTAppnyu=zeros(size(Nlist));
timeM=zeros(size(Nlist));
errNUFFT = zeros(size(Nlist));
errNUFFT1 = zeros(size(Nlist));
tol=1e-12;
num = 50;
iflag=-1;

for it = 1:length(Nlist)
    
    N = Nlist(it);
    
    x=sort(N*rand(N,1));
    k = (N-1)*rand(N,1);
    x1=2*pi*x/N;
    
    tic;
   
    nufftfun = nufftIII(k,x,iflag,N,100,tol);
    
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
    fk1 = nufft1d3(N,x1,c,iflag,tol,N,k);
    end
    timNUFFTAppnyu(it)=toc/num;
    
    
    [fhatM,ffun] = DeCom_NUFFT1D_III(c,x/N,k,tol);
    tic;
    for cnt = 1:num
        fhatM = ffun(c);
    end
    timeM(it) = toc/num;
    
    
    errNUFFT(it)=norm(fk1-fhatM,2)/norm(fk1,2);
    errNUFFT1(it)=norm(fk1-nufftc,2)/norm(fk1,2);
end

timNUFFTApp
%timNUFFTAppold
timNUFFTAppnyu
timeM
errNUFFT,
errNUFFT1
figure
loglog(Nlist,timNUFFTAppnyu);
hold on;
loglog(Nlist,timeM);
hold on;
loglog(Nlist,timNUFFTApp)
xlabel('N'),ylabel('time : s'),title('Comparison,1D type III,tol =1e-12'),legend('timeNYU','timeHZ','timeYZ','Location','northwest')
%errNUFFT1
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