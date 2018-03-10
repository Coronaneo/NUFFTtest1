function nufft3fun = nufft3II(x,iflag,ns,rt,tol)
% NUFFT3II: Computation of nonuniform FFT3 in R^2 - Type II.
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
    ns = ceil(max(max(x)))+1;%%±¾À´ÊÇmaxel
end

if nargin < 4
    rt = 175;
end

if nargin < 5
    tol = 1e-12;
end

n = size(x,1);

[k1,k2,k3] = ndgrid(0:ns-1,0:ns-1,0:ns-1);
k = [k1(:) k2(:) k3(:)];

fftconst = iflag*1i/ns*2*pi;

ratiofun = @(x,k)exp(fftconst*(x-round(x))*k');
[U,V] = lowrank(x,k,ratiofun,tol,rt,rt);

xsub = mod(round(x),ns)+1;
xxsub = sub2ind([ns ns ns],xsub(:,1),xsub(:,2),xsub(:,3));

nufft3fun = @(c)nufft3IIfun(c);

    function fft3c = nufft3IIfun(c)
        ncol = size(c,2);
        r = size(V,2);

        c = reshape( repmat(conj(V),1,ncol).*reshape(repmat(c,r,1), ns^3, r*ncol), ns,ns,ns,r*ncol);
        if iflag < 0
            fft3c = fft3(c);
        else
            fft3c = ifft3(c);
        end
        fft3c = reshape(fft3c,ns^3,r*ncol);
        fft3c = fft3c(xxsub,:);
        fft3c = squeeze(sum(reshape(repmat(U,1,ncol).*fft3c,n,r,ncol),2));
    end

end