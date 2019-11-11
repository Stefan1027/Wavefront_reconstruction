%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function, called NVWCER4SID, (which stands for Noise Variance Weighted
% Complex Exponential Reconstructor for Shearing Interferometer Data), which
% implements a version of the Itek exponential reconstructor---a version
% based on multiplication of the exponentials of i times the phase
% differences, i.e. on multiplication of what are called differential phasors
% (or DP's). It produces a set of complex values (called phasors)
% representing the reconstructed field, and associated variance like
% quantities. 
% 
% This function is a slight variation on the function called
% ComplexExpReconstructor, that was presented in TN-092. The variation
% consists in allowing the array representing the variance like quantity
% characteristic of the reconstructed wave function to be an output along
% with the array representing the reconstructed wave function. (The original
% ComplexExpReconstructor function only provided the array representing the
% reconstructed wave function as an output.)
% 
% This function assumes a set of (2^N+1)-by-(2^N) x-oriented DP's and a set
% of (2^N)-by-(2^N+1) y-oriented DP's, and produces a set of
% (2^N+1)-by-(2^N+1) phasors representing the reconstructed (complex) field.
% The DP's are input as unit magnitude quantities. Associated with each DP is
% a variance-like quantity (a "VLQ").
% 
% The reconstruction process consists of three parts, which are refered to
% here as the "reduction" part, the "simple-solve" part, and the
% "construction" part. The reduction phase is an iterative process, in each
% iteration of which a (nominally) factor-of-two smaller array of DP's are
% developed from the set of DP's developed in the preceeding iteration (or
% for the first iteration the DP's available as the inputs to this function).
% At the end of each iteration the thus far reduced set of DP's is saved
% before starting on the next iteration. In the N-th  iterations a pair of
% sets of DP's (along with associated VLQ's) of size 2-by-1 for the
% x-oriented DP's and of size 1-by-2 for the y-oriented DP's is produced.
% This is the last iteration of the reduction part of the process.  
% 
% In the  simple-solve part of the process a set of four phasors are
% calculated for the four corners of a square pattern, with the values of
% these four phasors being such that their differences (actually the product
% of one times the conjugate of the other) of the phasors for an adjacent
% pair of points of the square, represents a good approximation, in a least
% mean square error sense, to the DP values produced by the reduction part of
% this operation. This provides a 2-by-2 [or (2^n+1)-by-(2^n+1) for n=0] size
% array of phasor values. In addition to solving for the four phasors in the
% simple-solve part of the operation values of VLQ's are developed for each
% of these four phasors. The process is actually carried out with the DP's
% converted to phase differences and the with the results developed in terms
% of phases. (A value of 2*pi may be added or subtracted from one of the
% phase differences so as to minimize the closure error for the four phase
% differences.) At the end of the simple-solve the complex phasors for these
% phases are calculated.
% 
% In the construction part of this operation, in an iterative sequence of
% steps, a series of increasingly larger arrays of phasor values are
% developed. The first iteration starts with the 2-by-2 array developed by
% the simple-solve part of the operation and developes an array of phasors of
% size 3-by-3 [i.e. of size (2^n+1)-by-(2^n+1), for n=1]. In successive
% iterations the size of the array of phasors produced is  of size
% (2^n+1)-by-(2^n+1) for n=2, then for n=3, ... and eventually for n=N
% ---which is the size array of the desired final result. In the n-th
% iteration the phasor array at the start has a size of
% (2^(n-1)+1)-by-(2^(n-1)+1) and by the end of that iteration a phasor array
% of size (2^n+1)-by-(2^n+1) has been develped. The construction part of the
% operation is based on combining the set of phasor values available at the
% start of the iteration with the just higher spatial resolution set of DP's
% developed in the reduction part of the opteration. 
% 
% In the reduction part of the operation DP's are combined in two ways,
% in-series and in-parallel. The in-series combining occurs where a series of
% DP's  define a path from one point to an other point and we wish to develop
% a value for the DP between those two points. The in-parallel combining is
% where we wish to deal with the fact that there are several paths between
% two points, for each of which paths we have a DP value (the different DP
% values differing amongst themselves) and we wish to establish a
% reconciled/average/nominal value for the DP between those two points. The
% in-series DP value is obtained by simply multiplying the (complex) values
% of the DP's that define the path between the two points. Associated with
% this DP value will be a VLQ which is calculated as the sum of the VLQ's
% that are associated with each of the DP's on the path. The in-parallel DP
% value is calculated as sum of the (complex) values of the DP's for each of
% the several paths, each such DP value being multiplied by the inverse of
% the VLQ associated with that DP value before the summation, and then
% dividing that sum by the sum of the inverses of each of these VLQ's. The
% VLQ value associated with this in-parallel DP value result is calculated as
% the inverse of the sum of the inverses of the VLQ's associated with the DP
% values used in the in-parallel DP calculation.
% 
% In the simple-solve part of the process a least mean square error solution
% is found for the four phasor values which are the least incompatible with
% the four DP values. The concept of least incompatible takes account of the
% VLQ's associated with each of the four DP's (the calculations being
% conducted in terms of phases and phase differences, and not directly in
% terms of the complex phasors and DP's).
% 
% In each iterative cycle of the construction part of the operation a set of
% additional phasor values are to be developed from the previous set of
% phasors values combined with the relevant set of DP values. This process
% will involve both in-series and in-parellel sort of processes. The
% in-series process will differ from that in the reduction process only in
% that it will start from a phasor, join that with DP values, and produce a
% phasor value. As in the reduction part of the process the calcuation will
% simply require the multiplication of the relevant complex quantities. The
% VLQ associated with the resultant phasor will be simply the sum of the
% VLQ's associated with each of the complex quantities in the series. The
% in-parallel process will actually be simply the calculation of the weighted
% average of several calculated phasor values, all for the same point, to
% arrive at a reconciled/average/nominal value for the phasor at that point.
% This in-parallel phasor value is calculated as sum of the (complex) values
% of the several phasors, each such phasor value being multiplied (before the
% summation) by the inverse of the VLQ associated with that phasor value, and
% then dividing that sum by the sum of the inverses of each of these VLQ's.
% The VLQ value associated with this in-parallel phasor value result is
% calculated as the inverse of the sum of the inverses of the VLQ's
% associated with the phasor values used in the in-parallel phasor
% calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%          INPUT
%   DPx    (2^N+1,2^N)   = Differential phasor for displacement between
%                          array points in the x-direction, i.e. from one 
%                          column to the next column.
%   DPy    (2^N,2^N+1)   = Differential phasor for displacement between 
%                          array points in the y-direction, i.e. from one
%                          row to the next row. 
%   VLQx   (2^N+1,2^N)   = Quantity proportional to the nominal variance of
%                          the elements of PDc. This quantity has the same
%                          constant of proportionality as does Varr.
%   VLQy   (2^N,2^N+1)   = Quantity proportional to the nominal variance of
%                          the elements of PDr. This quantity has the same
%                          constant of proportionality as does Varc.  
% 
%          OUTPUTS
%   Phasor (2^N+1,2^N+1) = Reconstructed wave function. Phasor is a complex
%                          quantity.
%   VLQp   (2^N+1,2^N+1) = Variance like quantity characteristic of the
%                          variance of the quantity called Phasor. VLQp is
%                          a real quantity. 
%
%       [Phasor,VLQp]=NVWCER4SID(DPx,DPy,VLQx,VLQy)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phasor,VLQp]=NVWCER4SID(DPx,DPy,VLQx,VLQy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we check that the dimensions of the inputs are appropriate. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b]=size(DPx);
if a~=b+1
  error('DPx input is not (M+1)-by-M.')
end
N=log2(b);
if b~=2^round(N)
  error('DPx size is not (2^N+1)-by-(2^N).')
end
[c,d]=size(VLQx);
if c~=a | d~=b
  error('Size of VLQx does not match the size of DPx.')
end

[a,b]=size(DPy);
if a+1~=b
  error('DPy input is not M-by-(M+1).')
end
NN=log2(a);
if NN~=N
  error('DPy size is not (2^N)-by-(2^N+1).')
end
[c,d]=size(VLQy);
if c~=a | d~=b
  error('Size of VLQy does not match the size of DPy.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting the calculations we replace all infinite VLQ values by %
% a number that is 1e25 times larger than the largest finite VLQ value.  %
% This is done so as to avoid divide by zero operations in some of the   %
% computations. We also replace all zero DP values by a number that is   %
% 1e-10 smaller than the smallest non zero DP value.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii=find(isfinite(VLQx));
M=max(VLQx(ii));
ii=find(isinf(VLQx));
VLQx(ii)=1e10*M;
ii=find(isfinite(VLQy));
M=max(VLQy(ii));
ii=find(isinf(VLQy));
VLQy(ii)=1e10*M;

ax=abs(DPx);
ii=find(ax~=0);
M=min(ax(ii));
ii=find(ax==0);
DPx(ii)=1e-10*M;
ay=abs(DPy);
ii=find(ay~=0);
M=min(ay(ii));
ii=find(ay==0);
DPy(ii)=1e-10*M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We here start the REDUCTION process. We generate a set of DP's and         %
% associated VLQ's for reduced resolutions, denoted DPx##, DPy##, VLQx##,    %
% and VLQy## ---where ## is the power of two by which the size of the arrays %
% has been reduced relative to the input arrays. This is conducted in an     %
% iterative manner, reducing the size by successive factors of two till DPx  %
% and DPy are 2-by-1 and 1-by-2 respectively.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we initialize, setting the values for DPx##, DPy##, VLQx##, and %
% VLQy## (where ## denotes a two-digit representation of N) to the      %
% corresponding input arrays.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N<10
  eval(['DPx0' num2str(N) '=DPx;'])
  eval(['VLQx0' num2str(N) '=VLQx;'])
  eval(['DPy0' num2str(N) '=DPy;'])
  eval(['VLQy0' num2str(N) '=VLQy;'])
else
  eval(['DPx' num2str(N) '=DPx;'])
  eval(['VLQx' num2str(N) '=VLQx;'])
  eval(['DPy' num2str(N) '=DPy;'])
  eval(['VLQy' num2str(N) '=VLQy;'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In each cycle of the REDUCTION process we will start with a set of      %
% arrays, DPx***, DPy***, VLQx***, and VLQy***, where *** represents      % 
% "Old", and will generate a corresponding set of arrays of half their    %
% size for which *** represents "New". At the end of each cycle the "New" %
% values will be relabeled "Old" so as to prepare for start of the next   %
% cycle. First we set up the intial Old values.                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DPxOld=DPx;
VLQxOld=VLQx;
DPyOld=DPy;
VLQyOld=VLQy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we start the REDUCTION part of the exponential reconstruction process. %
% The index nu, or rather the value of N-nu, indicates the power of two      %
% corresponding to the size of the reduced arrays being developed in each    %
% cycle of the nu-loop.                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nu=1:N
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First we calculate the reduced resolution values for DPx and VLQx.   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  DPxNew=NaN*ones(2^(N-nu)+1,2^(N-nu));
  VLQxNew=NaN*ones(2^(N-nu)+1,2^(N-nu));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the first row.                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1
    for m=1:2^(N-nu)
      V1=VLQyOld(2*n-1,2*m-1)+VLQxOld(2*n,2*m-1) ...
                +VLQxOld(2*n,2*m)+VLQyOld(2*n-1,2*m+1);
      V2=VLQxOld(2*n-1,2*m-1)+VLQxOld(2*n-1,2*m);
      Vinv=1/V1+1/V2;
      VLQxNew(n,m)=1/Vinv;
      X1=DPyOld(2*n-1,2*m-1)*DPxOld(2*n,2*m-1)*DPxOld(2*n,2*m) ...
                *conj(DPyOld(2*n-1,2*m+1));
      X2=DPxOld(2*n-1,2*m-1)*DPxOld(2*n-1,2*m);
      DPxNew(n,m)=(X1/V1+X2/V2)/Vinv;
      if ~isfinite(DPxNew(n,m)); DPxNew(n,m)=1; end
      DPxNew(n,m)=DPxNew(n,m)/abs(DPxNew(n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For all but the first and last row.  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=2:2^(N-nu)
    for m=1:2^(N-nu)
      V1=VLQyOld(2*n-1,2*m-1)+VLQxOld(2*n,2*m-1) ...
                +VLQxOld(2*n,2*m)+VLQyOld(2*n-1,2*m+1);
      V2=VLQxOld(2*n-1,2*m-1)+VLQxOld(2*n-1,2*m);
      V3=VLQyOld(2*n-2,2*m-1)+VLQxOld(2*n-2,2*m-1) ...
                +VLQxOld(2*n-2,2*m)+VLQyOld(2*n-2,2*m+1);
      Vinv=1/V1+1/V2+1/V3;
      VLQxNew(n,m)=1/Vinv;
      X1=DPyOld(2*n-1,2*m-1)*DPxOld(2*n,2*m-1)*DPxOld(2*n,2*m) ...
                *conj(DPyOld(2*n-1,2*m+1));      
      X2=DPxOld(2*n-1,2*m-1)*DPxOld(2*n-1,2*m);
      X3=conj(DPyOld(2*n-2,2*m-1))*DPxOld(2*n-2,2*m-1)*DPxOld(2*n-2,2*m) ...
         *DPyOld(2*n-2,2*m+1);
      DPxNew(n,m)=(X1/V1+X2/V2+X3/V3)/Vinv;
      if ~isfinite(DPxNew(n,m)); DPxNew(n,m)=1; end
      DPxNew(n,m)=DPxNew(n,m)/abs(DPxNew(n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the last row.                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=2^(N-nu)+1
    for m=1:2^(N-nu)
      V2=VLQxOld(2*n-1,2*m-1)+VLQxOld(2*n-1,2*m);
      V3=VLQyOld(2*n-2,2*m-1)+VLQxOld(2*n-2,2*m-1) ...
                +VLQxOld(2*n-2,2*m)+VLQyOld(2*n-2,2*m+1);
      Vinv=1/V2+1/V3;
      VLQxNew(n,m)=1/Vinv;
      X2=DPxOld(2*n-1,2*m-1)*DPxOld(2*n-1,2*m);
      X3=conj(DPyOld(2*n-2,2*m-1))*DPxOld(2*n-2,2*m-1)*DPxOld(2*n-2,2*m) ...
         *DPyOld(2*n-2,2*m+1);
      DPxNew(n,m)=(X2/V2+X3/V3)/Vinv;
      if ~isfinite(DPxNew(n,m)); DPxNew(n,m)=1; end
      DPxNew(n,m)=DPxNew(n,m)/abs(DPxNew(n,m));
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we calculate the reduced resolution values for DPy and VLQy.   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  DPyNew=NaN*ones(2^(N-nu),2^(N-nu)+1);
  VLQyNew=NaN*ones(2^(N-nu),2^(N-nu)+1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the first column.               %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for m=1
    for n=1:2^(N-nu)
      V1=VLQxOld(2*n-1,2*m-1)+VLQyOld(2*n-1,2*m) ...
                +VLQyOld(2*n,2*m)+VLQxOld(2*n+1,2*m-1);
      V2=VLQyOld(2*n-1,2*m-1)+VLQyOld(2*n,2*m-1);
      Vinv=1/V1+1/V2;
      VLQyNew(n,m)=1/Vinv;
      X1=DPxOld(2*n-1,2*m-1)*DPyOld(2*n-1,2*m)*DPyOld(2*n,2*m) ...
                *conj(DPxOld(2*n+1,2*m-1));
      X2=DPyOld(2*n-1,2*m-1)*DPyOld(2*n,2*m-1);
      DPyNew(n,m)=(X1/V1+X2/V2)/Vinv;
      if ~isfinite(DPyNew(n,m)); DPyNew(n,m)=1; end
      DPyNew(n,m)=DPyNew(n,m)/abs(DPyNew(n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For all but the first and last columns. %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for m=2:2^(N-nu)
    for n=1:2^(N-nu)
      V1=VLQxOld(2*n-1,2*m-1)+VLQyOld(2*n-1,2*m) ...
                +VLQyOld(2*n,2*m)+VLQxOld(2*n+1,2*m-1);
      V2=VLQyOld(2*n-1,2*m-1)+VLQyOld(2*n,2*m-1);
      V3=VLQxOld(2*n-1,2*m-2)+VLQyOld(2*n-1,2*m-2) ...
                +VLQyOld(2*n,2*m-2)+VLQxOld(2*n+1,2*m-2);
      Vinv=1/V1+1/V2+1/V3;
      VLQyNew(n,m)=1/Vinv;
      X1=DPxOld(2*n-1,2*m-1)*DPyOld(2*n-1,2*m)*DPyOld(2*n,2*m) ...
                *conj(DPxOld(2*n+1,2*m-1));
      X2=DPyOld(2*n-1,2*m-1)*DPyOld(2*n,2*m-1);
      X3=conj(DPxOld(2*n-1,2*m-2))*DPyOld(2*n-1,2*m-2)*DPyOld(2*n,2*m-2) ...
         *DPxOld(2*n+1,2*m-2);
      DPyNew(n,m)=(X1/V1+X2/V2+X3/V3)/Vinv;
      if ~isfinite(DPyNew(n,m)); DPyNew(n,m)=1; end
      DPyNew(n,m)=DPyNew(n,m)/abs(DPyNew(n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the last column.                %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for m=2^(N-nu)+1
    for n=1:2^(N-nu)
      V2=VLQyOld(2*n-1,2*m-1)+VLQyOld(2*n,2*m-1);
      V3=VLQxOld(2*n-1,2*m-2)+VLQyOld(2*n-1,2*m-2) ...
                +VLQyOld(2*n,2*m-2)+VLQxOld(2*n+1,2*m-2);
      Vinv=1/V2+1/V3;
      VLQyNew(n,m)=1/Vinv;
      X2=DPyOld(2*n-1,2*m-1)*DPyOld(2*n,2*m-1);
      X3=conj(DPxOld(2*n-1,2*m-2))*DPyOld(2*n-1,2*m-2)*DPyOld(2*n,2*m-2) ...
         *DPxOld(2*n+1,2*m-2);
      DPyNew(n,m)=(X2/V2+X3/V3)/Vinv;
      if ~isfinite(DPyNew(n,m)); DPyNew(n,m)=1; end
      DPyNew(n,m)=DPyNew(n,m)/abs(DPyNew(n,m));
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Before we can start the next cycle of the nu-loop we have to relabel %
  % our New-values, now calling them Old-values --- for use as the       %
  % Old-values in the next cycle, and also have to store the New values  %
  % in numbered arrays for use latter in the CONSTRUCTION part of the    %
  % exponential reconstruction process.                                  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  VLQxOld=VLQxNew;
  VLQyOld=VLQyNew;
  DPxOld=DPxNew;
  DPyOld=DPyNew;

  if N-nu<10
    eval(['DPx0' num2str(N-nu) '=DPxNew;'])
    eval(['VLQx0' num2str(N-nu) '=VLQxNew;'])
    eval(['DPy0' num2str(N-nu) '=DPyNew;'])
    eval(['VLQy0' num2str(N-nu) '=VLQyNew;'])
  else
    eval(['DPx' num2str(N-nu) '=DPxNew;'])
    eval(['VLQx' num2str(N-nu) '=VLQxNew;'])
    eval(['DPy' num2str(N-nu) '=DPyNew;'])
    eval(['VLQy' num2str(N-nu) '=VLQyNew;'])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have now completed the generation of all the reduced resolution DP's  %
% and associated VLQ's. The REDUCTION part of exponential reconstruction   %
% process is now complete and we are ready to start the SIMPLE-SOLVE part  %
% of the process. This will be conducted as a weighted least mean square   %
% error PHASE reconstruction, with the PHASOR being calculted from the     %
% PHASE.                                                                   %
%                                                                          %
% We have reduced the DP's and associated VLQ's to a pair of 2-by-1 arrays %
% and a pair of 1-by-2 arrays. We first convert the DP's (differential     %
% phasors) to phase differences. With this we start the calculation of the %
% implied 2-by-2 array of phases and associated variances, and then of the %
% corresponding complex phasors.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta=angle([DPx00.' DPy00])';
curl=Delta(1)-Delta(2)+Delta(3)-Delta(4);
if curl>pi
  Delta(1)=Delta(1)-2*pi;
end
if curl<-pi
  Delta(1)=Delta(1)+2*pi;
end

R=diag(1./[VLQx00' VLQy00]);
Gamma=[-1 1 0; -1 -1 -2; -1 0 1; -1 -2 -1];
H=inv(Gamma'*R*Gamma)*(Gamma'*R);
mu=[];
for k=1:4
  mu=[mu -H(1,k)-H(2,k)-H(3,k)];
end
H=[H; mu];
phi=reshape(H*Delta,2,2)';
PhasorOld=exp(i*phi);
VLQpOld=reshape(diag(H*inv(R)*H'),2,2)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This completes the solution for the phases of the 2-by-2 reduced size    %
% arrays. We are now ready to start the iterative part of the CONSTRUCTION %
% portion of the  exponential reconstruction process, building up the      %
% size of the array of phasors and associated variances in factors of      %
% two at each cycle of thenu-loop.                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nu=1:N
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % At the start of each cycle of the nu-loop we have the phasor, denoted  % 
  % by PhasorOld, and its variance, denoted by VLQpOld, defined on         %
  % (2^(nu-1)+1)-by-(2^(nu-1)+1) size arrays. In the cycle we will develop %
  % new values for the phasor and its variance, denoted by PhasorNew and   %
  % VLQpNew, on (2^nu+1)-by-(2^nu+1) size arrays. To accomplish this BUILD %
  % portion of the exponential reconstructor, improving the resolution we  %
  % will make use of the higher resolution phase difference data stored in %
  % DPx##, DPy##, VLQx##, and VLQy## (where ## is the two digit            %
  % representation of the value of nu), calling the recovered values       %
  % DPxOld, DPyOld, VLQxOld, and VLQyOld respectively. Just before         %
  % completing each cycle of the nu-loop we will relabel the PhasorNew and %
  % VLQpNew quantities as PhasorOld and VLQpOld, so as to be ready to      %
  % start the new cycle.                                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First we recall the higher resolution differential phasor data. %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nu<10
    eval(['DPxOld=DPx0' num2str(nu) ';'])
    eval(['DPyOld=DPy0' num2str(nu) ';'])
    eval(['VLQxOld=VLQx0' num2str(nu) ';'])
    eval(['VLQyOld=VLQy0' num2str(nu) ';'])
  else
    eval(['DPxOld=DPx' num2str(nu) ';'])
    eval(['DPyOld=DPy' num2str(nu) ';'])
    eval(['VLQxOld=VLQx' num2str(nu) ';'])
    eval(['VLQyOld=VLQy' num2str(nu) ';'])
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we start the development of PhasorNew and VLQpNew. %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  PhasorNew=NaN*ones(2^nu+1);
  VLQpNew=NaN*ones(2^nu+1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % First we insert into these two arrays the values that we know from % 
  % before the start of this cycle of the nu-loop.                     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)+1
    for m=1:2^(nu-1)+1
      PhasorNew(2*n-1,2*m-1)=PhasorOld(n,m);
      VLQpNew(2*n-1,2*m-1)=VLQpOld(n,m);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % The phasor values we have just filled into the PhasorNew array  %
  % can be considered to form a set of (2^(nu-1))-by-(2^(nu-1))     %
  % boxes, each with nine points of which only the four corners are %
  % filled in. We now proceedto calculate/fill-in the values for    %
  % PhasorNew and VLQpNew for the CENTERS of each of these boxes.   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)
    for m=1:2^(nu-1)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Starting from the lower left corner of the box. %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      v1=VLQxOld(2*n-1,2*m-1)+VLQyOld(2*n-1,2*m);
      v2=VLQyOld(2*n-1,2*m-1)+VLQxOld(2*n,2*m-1);
      x1=DPxOld(2*n-1,2*m-1)*DPyOld(2*n-1,2*m);
      x2=DPyOld(2*n-1,2*m-1)*DPxOld(2*n,2*m-1);
      v=1/(1/v1+1/v2);
      x=(x1/v1+x2/v2)*v;
      if ~isfinite(x); x=1; end
      x=x/abs(x);
      V1=v+VLQpNew(2*n-1,2*m-1);
      X1=x*PhasorNew(2*n-1,2*m-1);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Starting from the lower right corner of the box. %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      v1=VLQxOld(2*n-1,2*m)+VLQyOld(2*n-1,2*m);
      v2=VLQyOld(2*n-1,2*m+1)+VLQxOld(2*n,2*m);
      x1=conj(DPxOld(2*n-1,2*m))*DPyOld(2*n-1,2*m);
      x2=DPyOld(2*n-1,2*m+1)*conj(DPxOld(2*n,2*m));
      v=1/(1/v1+1/v2);
      x=(x1/v1+x2/v2)*v;
      if ~isfinite(x); x=1; end
      x=x/abs(x);
      V2=v+VLQpNew(2*n-1,2*m+1);
      X2=x*PhasorNew(2*n-1,2*m+1);

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Starting from the upper left corner of the box. %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      v1=VLQxOld(2*n+1,2*m-1)+VLQyOld(2*n,2*m);
      v2=VLQyOld(2*n,2*m-1)+VLQxOld(2*n,2*m-1);
      x1=DPxOld(2*n+1,2*m-1)*conj(DPyOld(2*n,2*m));
      x2=conj(DPyOld(2*n,2*m-1))*DPxOld(2*n,2*m-1);
      v=1/(1/v1+1/v2);
      x=(x1/v1+x2/v2)*v;
      if ~isfinite(x); x=1; end
      x=x/abs(x);
      V3=v+VLQpNew(2*n+1,2*m-1);
      X3=x*PhasorNew(2*n+1,2*m-1);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Starting from the upper right corner of the box. %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      v1=VLQxOld(2*n+1,2*m)+VLQyOld(2*n,2*m);
      v2=VLQyOld(2*n,2*m+1)+VLQxOld(2*n,2*m);
      x1=conj(DPxOld(2*n+1,2*m))*conj(DPyOld(2*n,2*m));
      x2=conj(DPyOld(2*n,2*m+1))*conj(DPxOld(2*n,2*m));
      v=1/(1/v1+1/v2);
      x=(x1/v1+x2/v2)*v;
      if ~isfinite(x); x=1; end
      x=x/abs(x);
      V4=v+VLQpNew(2*n+1,2*m+1);
      X4=x*PhasorNew(2*n+1,2*m+1);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Now combining the contributions from the four corners of the   %
      % box (via the two paths from each corner to the center of the   %
      % box) to get the phasor, and its variance, at the center of the %
      % box.                                                           %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      VLQpNew(2*n,2*m)=1/(1/V1+1/V2+1/V3+1/V4);
      PhasorNew(2*n,2*m)=(X1/V1+X2/V2+X3/V3+X4/V4)*VLQpNew(2*n,2*m);
      if ~isfinite(PhasorNew(2*n,2*m)); PhasorNew(2*n,2*m)=1; end
      if PhasorNew(2*n,2*m)==0
        disp('zero value for PhasorNew')
        keyboard
      end
      PhasorNew(2*n,2*m)=PhasorNew(2*n,2*m)/abs(PhasorNew(2*n,2*m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This complets the calculation of the phasor and its variance at the %
  % CENTER of each of the (2^(nu-1))-by-(2^(nu-1)) boxes. We now start  %
  % calculation of the phasor and its variance for the position at the  %
  % center of the BOX EDGES on the top and on the right for all the     %
  % boxes except the top-most and the right-most boxes. (This will take % 
  % care of the positions at the center of each of the box edges except %
  % those for edges which are part of the outside edges of the large    %
  % box enclossing the entire array.)                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the upper edges of each box     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)-1
    for m=1:2^(nu-1)
      V1=VLQpNew(2*n,2*m)+VLQyOld(2*n,2*m);              % Start below
      X1=PhasorNew(2*n,2*m)*DPyOld(2*n,2*m);
      X1=X1/abs(X1);
      V2=VLQpNew(2*n+1,2*m-1)+VLQxOld(2*n+1,2*m-1);      % Start left
      X2=PhasorNew(2*n+1,2*m-1)*DPxOld(2*n+1,2*m-1);
      X2=X2/abs(X2);
      V3=VLQpNew(2*n+2,2*m)+VLQyOld(2*n+1,2*m);          % Start above
      X3=PhasorNew(2*n+2,2*m)*conj(DPyOld(2*n+1,2*m));
      X3=X3/abs(X3);
      V4=VLQpNew(2*n+1,2*m+1)+VLQxOld(2*n+1,2*m);        % Start right
      X4=PhasorNew(2*n+1,2*m+1)*conj(DPxOld(2*n+1,2*m));
      X4=X4/abs(X4);
      VLQpNew(2*n+1,2*m)=1/(1/V1+1/V2+1/V3+1/V4);
      PhasorNew(2*n+1,2*m)=(X1/V1+X2/V2+X3/V3+X4/V4)*VLQpNew(2*n+1,2*m);
      if ~isfinite(PhasorNew(2*n+1,2*m)); PhasorNew(2*n+1,2*m)=1; end
      PhasorNew(2*n+1,2*m)=PhasorNew(2*n+1,2*m)/abs(PhasorNew(2*n+1,2*m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the right edges of each box     %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)
    for m=1:2^(nu-1)-1
      V1=VLQpNew(2*n-1,2*m+1)+VLQyOld(2*n-1,2*m+1);      % Start below
      X1=PhasorNew(2*n-1,2*m+1)*DPyOld(2*n-1,2*m+1);
      V2=VLQpNew(2*n,2*m)+VLQxOld(2*n,2*m);              % Start left
      X2=PhasorNew(2*n,2*m)*DPxOld(2*n,2*m);
      V3=VLQpNew(2*n+1,2*m+1)+VLQyOld(2*n,2*m+1);        % Start above
      X3=PhasorNew(2*n+1,2*m+1)*conj(DPyOld(2*n,2*m+1));
      V4=VLQpNew(2*n,2*m+2)+VLQxOld(2*n,2*m+1);          % Start right
      X4=PhasorNew(2*n,2*m+2)*conj(DPxOld(2*n,2*m+1));
      VLQpNew(2*n,2*m+1)=1/(1/V1+1/V2+1/V3+1/V4);
      PhasorNew(2*n,2*m+1)=(X1/V1+X2/V2+X3/V3+X4/V4)*VLQpNew(2*n,2*m+1);
      if ~isfinite(PhasorNew(2*n,2*m+1)); PhasorNew(2*n,2*m+1)=1; end
      PhasorNew(2*n,2*m+1)=PhasorNew(2*n,2*m+1)/abs(PhasorNew(2*n,2*m+1));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This completes the evaluation of the phasor and its variance at the %
  % center of the EDGES of all the elemental boxes escept for the edges %
  % that bound the array. We now start the evaluation for those four    %
  % sets of elemental BOUNDING EDGES.                                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the bottom edge of the array %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1
    for m=1:2^(nu-1)
      V1=VLQpNew(n,2*m-1)+VLQxOld(n,2*m-1);           % Start left
      X1=PhasorNew(n,2*m-1)*DPxOld(n,2*m-1);
      V2=VLQpNew(n+1,2*m)+VLQyOld(n,2*m);             % Start above
      X2=PhasorNew(n+1,2*m)*conj(DPyOld(n,2*m));
      V3=VLQpNew(n,2*m+1)+VLQxOld(n,2*m);             % Start right
      X3=PhasorNew(n,2*m+1)*conj(DPxOld(n,2*m));
      VLQpNew(n,2*m)=1/(1/V1+1/V2+1/V3);
      PhasorNew(n,2*m)=(X1/V1+X2/V2+X3/V3)*VLQpNew(n,2*m);
      if ~isfinite(PhasorNew(n,2*m)); PhasorNew(n,2*m)=1; end
      PhasorNew(n,2*m)=PhasorNew(n,2*m)/abs(PhasorNew(n,2*m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the upper edge of the array %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=2^nu+1
    for m=1:2^(nu-1)
      V1=VLQpNew(n,2*m+1)+VLQxOld(n,2*m);             % Start right
      X1=PhasorNew(n,2*m+1)*conj(DPxOld(n,2*m));
      V2=VLQpNew(n-1,2*m)+VLQyOld(n-1,2*m);           % Start below
      X2=PhasorNew(n-1,2*m)*DPyOld(n-1,2*m);
      V3=VLQpNew(n,2*m-1)+VLQx(n,2*m-1);              % Start left
      X3=PhasorNew(n,2*m-1)*DPxOld(n,2*m-1);
      VLQpNew(n,2*m)=1/(1/V1+1/V2+1/V3);
      PhasorNew(n,2*m)=(X1/V1+X2/V2+X3/V3)*VLQpNew(n,2*m);
      if ~isfinite(PhasorNew(n,2*m)); PhasorNew(n,2*m)=1; end
      PhasorNew(n,2*m)=PhasorNew(n,2*m)/abs(PhasorNew(n,2*m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the left edge of the array %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)
    for m=1
      V1=VLQpNew(2*n+1,m)+VLQyOld(2*n,m);             % Start above
      X1=PhasorNew(2*n+1,m)*conj(DPyOld(2*n,m));
      V2=VLQpNew(2*n,m+1)+VLQxOld(2*n,m);             % Start right
      X2=PhasorNew(2*n,m+1)*conj(DPxOld(2*n,m));
      V3=VLQpNew(2*n-1,m)+VLQyOld(2*n-1,m);           % Start below
      X3=PhasorNew(2*n-1,m)*DPyOld(2*n-1,m);
      VLQpNew(2*n,m)=1/(1/V1+1/V2+1/V3);
      PhasorNew(2*n,m)=(X1/V1+X2/V2+X3/V3)*VLQpNew(2*n,m);
      if ~isfinite(PhasorNew(2*n,m)); PhasorNew(2*n,m)=1; end
      PhasorNew(2*n,m)=PhasorNew(2*n,m)/abs(PhasorNew(2*n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % For the right edge of the array %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for n=1:2^(nu-1)
    for m=2^nu+1
      V1=VLQpNew(2*n-1,m)+VLQyOld(2*n-1,m);           % Start below
      X1=PhasorNew(2*n-1,m)*DPyOld(2*n-1,m);
      V2=VLQpNew(2*n,m-1)+VLQxOld(2*n,m-1);           % Start left
      X2=PhasorNew(2*n,m-1)*DPxOld(2*n,m-1);
      V3=VLQpNew(2*n+1,m)+VLQyOld(2*n,m);             % Start above
      X3=PhasorNew(2*n+1,m)*conj(DPyOld(2*n,m));
      VLQpNew(2*n,m)=1/(1/V1+1/V2+1/V3);
      PhasorNew(2*n,m)=(X1/V1+X2/V2+X3/V3)*VLQpNew(2*n,m);
      if ~isfinite(PhasorNew(2*n,m)); PhasorNew(2*n,m)=1; end
      PhasorNew(2*n,m)=PhasorNew(2*n,m)/abs(PhasorNew(2*n,m));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This completes the generation of the phasor and associated variance    %
  % like quantity, PhasorNew and VLQpNew respectively, for the current     %
  % cycle of the nu-loop. The size of these arrays are                     %
  % (2^nu+1)-by-(2^nu+1). We now rename these variables as PhasorOld and   %
  % VLQpOld in preparation for the start of the next cycle of the nu-loop. %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  PhasorOld=PhasorNew;
  VLQpOld=VLQpNew;
end

Phasor=PhasorOld;
VLQp=VLQpOld;

% keyboard
