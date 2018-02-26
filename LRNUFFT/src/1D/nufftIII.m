function nufftfun = nufftIII(k,x,iflag,ns,rt,tol)
% NUFFTIII: Computation of nonuniform FFT in R^1 - Type III.
%   gfun = NUFFTIII(k,x,iflag,ns,rt,tol) provides the fast algorithm for
%   nufft of type II. k is the location of targets on interval [0,ns-1], x
%   is the location of sources on interval [0,ns-1], iflag determines the
%   sign of FFT, and ns is the number of Fourier modes computed. The
%   function returns a function that fast evaluate nufft.
%
%   NUFFT of type III is defined as follows.
%                 ns
%     g(k(i)) = w sum c(x(j)) exp(+/-1i 2pi * k(i) x(j) / n) 
%                 j=1
%   for 0 <= k(i) < ns and 0 <= x(j) < ns.
%   If (iflag .ge.0) the + sign is used in the exponential and a = 1/ns.
%   If (iflag .lt.0) the - sign is used in the exponential and a = 1.
%
%   See also FFT, NUFFTI, NUFFTII.


if nargin < 3
    iflag = -1;
else
    iflag = sign(iflag);
end

if nargin < 4
    ns = max(ceil(max(x)),ceil(max(k)))+1;
end

if nargin < 5
    rt = 100;
end

if nargin < 6
    tol = 1e-12;
end

fftconst = iflag*1i/ns*2*pi;

ratiofun = @(k,x)exp(fftconst*(k*x'-round(k)*round(x)'));
[U,V] = lowrank(k,x,ratiofun,tol,rt,rt);

xsub = mod(round(x),ns)+1;
spPerm = sparse(xsub,1:ns,ones(1,ns),ns,ns);
ksub = mod(round(k),ns)+1;

r = size(V,2);

nufftfun = @(c)nufftIIIfun(c);

    function fftc = nufftIIIfun(c)
        [n,ncol] = size(c);

        c = repmat(conj(V),1,ncol).*reshape(repmat(c,r,1),n,r*ncol);
        c = spPerm*c;
        if iflag < 0
            fftc = fft(c);
        else
            fftc = ifft(c);
        end
        fftc = fftc(ksub,:);
        fftc = squeeze(sum(reshape(repmat(U,1,ncol).*fftc,n,r,ncol),2));
    end

end