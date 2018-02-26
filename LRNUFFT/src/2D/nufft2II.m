function nufft2fun = nufft2II(x,iflag,ns,rt,tol)
% NUFFT2II: Computation of nonuniform FFT2 in R^2 - Type II.
%   gfun = NUFFT2II(x,iflag,ns,rt,tol) provides the fast algorithm for
%   nufft2 of type II. x is the location of sources on interval
%   [0,ns-1]x[0,ns-1], iflag determines the sign of FFT2, and ns is the
%   number of Fourier modes computed. The function returns a function that
%   fast evaluate nufft2.
%
%   NUFFT2 of type II is defined as follows.
%                            ns
%     f(x(j,:)) = w   sum   g(k1,k2) exp(+/-i 2pi * [k1 k2] x(j,:)' / ns) 
%                          k1,k2=1
%   for 1 <= k1,k2 <= ns and 1 <= x(j,1),x(j,2) <= ns.
%   If (iflag .ge.0) the + sign is used in the exponential and a = 1/ns^2.
%   If (iflag .lt.0) the - sign is used in the exponential and a = 1.
%
%   See also FFT, FFT2, NUFFTII, NUFFT2I, NUFFT2III.


if nargin < 2
    iflag = 1;
else
    iflag = sign(iflag);
end

if nargin < 3
    ns = ceil(maxel(x))+1;
end

if nargin < 4
    rt = 135;
end

if nargin < 5
    tol = 1e-12;
end

n = size(x,1);

[k1,k2] = ndgrid(0:ns-1);
k = [k1(:) k2(:)];

fftconst = iflag*1i/ns*2*pi;

ratiofun = @(x,k)exp(fftconst*(x-round(x))*k');
[U,V] = lowrank(x,k,ratiofun,tol,rt,rt);

xsub = mod(round(x),ns)+1;
xxsub = sub2ind([ns ns],xsub(:,1),xsub(:,2));

nufft2fun = @(c)nufft2IIfun(c);

    function fft2c = nufft2IIfun(c)
        ncol = size(c,2);
        r = size(V,2);

        c = reshape( repmat(conj(V),1,ncol) ...
            .*reshape(repmat(c,r,1), ns^2, r*ncol), ns, ns, r*ncol);
        if iflag < 0
            fft2c = fft2(c);
        else
            fft2c = ifft2(c);
        end
        fft2c = reshape(fft2c,ns^2,r*ncol);
        fft2c = fft2c(xxsub,:);
        fft2c = squeeze(sum(reshape(repmat(U,1,ncol).*fft2c,n,r,ncol),2));
    end

end