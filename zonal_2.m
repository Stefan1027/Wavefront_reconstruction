function [phases] = zonal_2(Fx,Fy)
h=1;
[a,b]=size(Fx);
[c,d]=size(Fy);

%Errors to determine if input data represents a square value or not
if c~=a | d~=a
    error('Size of Fx differs from Fy');
end
if a~=b
    error('Fx is not square');
end
if c~=d
    error('Fy is not square');
end

N=a;
S=[reshape(Fx',1,N^2) reshape(Fy',1,N^2)]';
A=formA(N);
D=formD(N);
phases=pinv(A)*D*S;
% [U,A2,V]=svd(A,0);
% A2=pinv(A2);
% phases=V*A2*U?*D*S;
phases=reshape(phases,N,N)'./h;
end

% Derived from (Dai 2008)
function A=formA(N)
A=zeros(2*N*(N-1),N^2);
for i=1:N
    for j=1:(N-1)
        A((i-1)*(N-1)+j,(i-1)*N+j)=-1;
        A((i-1)*(N-1)+j,(i-1)*N+j+1)=1;
        A((N+i-1)*(N-1)+j,i+(j-1)*N)=-1;
        A((N+i-1)*(N-1)+j,i+j*N)=1;
    end
end
end

function D=formD(N)
D=zeros(2*N*(N-1),2*N^2);
for i=1:N
    for j=1:(N-1)
        D((i-1)*(N-1)+j,(i-1)*N+j)=0.5;
        D((i-1)*(N-1)+j,(i-1)*N+j+1)=0.5;
        D((N+i-1)*(N-1)+j,N*(N+j-1)+i)=0.5;
        D((N+i-1)*(N-1)+j,N*(N+j)+i)=0.5;
    end
end
end
