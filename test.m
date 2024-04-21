x=[1,2,3];
n=length(x);%location of rxx(0);
m=n-1;% m+1<=n;
rx=xcorr(x);%length of rx is 2n-1;
Rxx=toeplitz(rx(n:n+m))/n