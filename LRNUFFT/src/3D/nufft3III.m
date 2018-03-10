function nufft3fun = nufft3III(k,x,iflag,ns,rt,tol)
% NUFFT2III: Computation of nonuniform FFT2 in R^2 - Type III.
%   gfun = NUFFT2III(k,x,iflag,ns,rt,tol) provides the fast algorithm for
%   nufft2 of type II. k is the location of targets on interval [0,ns-1]x[0,ns-1], x
%   is the location of sources on interval [0,ns-1]x[0,ns-1], iflag determines the
%   sign of FFT2, and ns is the number of Fourier modes computed. The
%   function returns a function that fast evaluate nufft2.
%
%   NUFFT2 of type III is defined as follows.
%                   n
%     g(k(i,:)) = w sum c(x(j,:)) exp(+/-1i 2pi * k(i,:) x(j,:)' / ns) 
%                   j=1
%   for 0 <= k(i,1),k(i,2) < ns and 0 <= x(j,1),x(j,2) < ns.
%   If (iflag .ge.0) the + sign is used in the exponential and a = 1/ns^2.
%   If (iflag .lt.0) the - sign is used in the exponential and a = 1.
%
%   See also FFT, FFT2, NUFFTIII, NUFFT2I, NUFFT2II.


if nargin < 3
    iflag = -1;
else
    iflag = sign(iflag);
end

if nargin < 4
    ns = max(ceil(max(max(x))),ceil(max(max(k))))+1;
end

if nargin < 5
    rt = 135;
end

if nargin < 6
    tol = 1e-12;
end

fftconst = iflag*1i/ns*2*pi;

ratiofun = @(k,x)exp(fftconst*(k-round(k))*x');
[U,V] = lowrank(k,x,ratiofun,tol,rt,rt);

ksub = mod(round(k),ns)+1;

nufft3Ifun = nufft3I(x,iflag,ns,rt,tol);
kksub = sub2ind([ns ns ns],ksub(:,1),ksub(:,2),ksub(:,3));

nufft3fun = @(c)nufft3IIIfun(c);

    function fft3c = nufft3IIIfun(c)
        [n,ncol] = size(c);
        r = size(V,2);
        
        c = repmat(conj(V),1,ncol).*reshape(repmat(c,r,1), n, r*ncol);
        fft3c = nufft3Ifun(c);
        fft3c = fft3c(kksub,:);
        fft3c = squeeze(sum(reshape(repmat(U,1,ncol).*fft3c,n,r,ncol),2));
    end

end