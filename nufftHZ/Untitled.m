format long
nj=65536;
ms=65536;
xj=(1:nj)'*pi/nj/4;
cj=exp(1i*(1:nj)/nj)';
k=-ms/2:(ms/2-1);
%fun = @(k,x) exp(-2*pi*i*k*(x-(mod(round(x*nj),nj)/nj))');
%K=13;
tol=1e-12;
%[U,V] = lowrank(k(:),xj(:),fun,tol,5*K,K);
%R=conj(V);
%r=size(V,2);
[fhat,ffun,L,R,Id,M1] = DeCom_NUFFT1D_I(cj,xj,k,tol);
%idx=mod(round(xj*nj),nj)+1;
%Id=sparse(idx,1:nj,ones(1,nj),nj,nj);
r=size(R,2)
fk1 = nufft1d1(nj,2*pi*xj,cj,-1,tol,nj)*nj;
norm(fhat-fk1)/norm(fhat)
M1r=real(M1);
M1i=real(-1i*M1);
%norm(M1r+1i*M1i-M1)
fid=fopen('M1rh.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M1r);
fclose(fid);
fid=fopen('M1ih.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M1i);
fclose(fid);



fhat1=fhat;
[fhat,ffun,L,R,Idx,M2] = DeCom_NUFFT1D_II(fhat1,xj,k,tol);
cj1 = nufft1d2(nj,2*pi*xj,-1,tol,nj,fhat1);
norm(cj1-fhat)/norm(fhat)
M2r=real(M2);
M2i=real(-i*M2);
%norm(M2r+i*M2i-M2)
fid=fopen('M2rh.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M2r);
fclose(fid);
fid=fopen('M2ih.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M2i);
fclose(fid);

[fhat,ffun,L,R,Idk,M3] = DeCom_NUFFT1D_III(cj,xj,k,tol);
fk1 = nufft1d3(nj,2*pi*xj,cj,-1,tol,nj,k);
norm(fhat-fk1)/norm(fhat)
M3r=real(M3);
M3i=real(-i*M3);
%norm(M3r+i*M3i-M3)
fid=fopen('M3rh.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M3r);
fclose(fid);
fid=fopen('M3ih.txt','w');
fprintf(fid,'%12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\r\n',M3i);
fclose(fid);

