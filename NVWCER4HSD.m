function [E,phi,VLQ]=NVWCER4HSD(dpx,dpy,sigsqx,sigsqy,mask,amp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function, called NVWCER4HSD, (which stands for Noise Variance Weighted
% Complex Exponential Reconstructor for Hartmann Sensor Data) which
% implements a version of the complex exponential reconstructor that will
% work on phase difference measurement values provided by a Hartmann wave
% front distortion sensor. The operation of this function is based on
% separating the phase difference values, reformed as sums and differences,
% into two sets (or sub spaces), in each of which the phase difference data
% is aranged in a Hutchin-geometry---so the NVWCER4SID function can be
% applied. This provides two reconstructed wave arrays of wave function data,
% one in each of the two sub spaces. Applying these results to the original
% sample space places the results as two interleaved sets of results, which
% have to be brought into phase with each other.
% 
% This function is organized to operate in three phases. In the first phase
% the Hartmann sensor data is arranged to support a Hutchin-geometry for the
% first sub space, and reconstruction of the complex wave function values for
% the first sub space is carried out. In the second phase the Hartmann sensor
% data is arranged to support a Hutchin-geometry for the second sub space,
% and reconstruction of the complex wave function values for the second sub
% space is carried out. In the third phase the reconstructed wave function
% results developed for the first sub space and the reconstructed wave
% function results developed for the second sub space are combined in an
% interleaved form, and the two have their phases adjusted so as to minimize
% any egg-crating type discrepancy.
% 
% In the first and second phases of the operation of this function-program,
% in order to develop the Hutchin-geometry array of phase differences needed
% as inputs to the NVWCER4SID function phase difference values are needed for
% positions such that the values can not be obtained from the available
% Hartmann sensor data. These "unobtainable" values are set to zero and have
% an infinite variance assigned to them. In the third phase of the operation
% of this function-program, in transfering the complex wave function values
% from the first and second sub spaces to the original sample space (on which
% the reconstructed wave function results are to be developed) some of the
% values in the first and second sub spaces will fall outside the range of
% interest of the original sample space; these values will be discarded.
% 
% As an option the estimated amplitude of the reconstructed wave funcion at
% each of the corners of the Hartmann sensor's sub apertures can be provided
% as an input. If this input is provided the amplitudes of the reconstructed
% wave function will be adjusted to match these values. If not provided then
% the reconstructed wave function will have unity amplitude at all points.
% 
% The input called mask indicates the positions corresponding to points that
% are within the aperture. The elements (positions) of mask correspond to
% corners of the Hartmann sensor's sub apertures, for each of the elements of
% which one of the elements of dpx, dpy, sigsqx,
% and sigsqy correspond. Any element of sigsqx and of sigsqy which has a
% value of infinity corresponds to a sub aperture for which the four corners
% all have mask element values of zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Adapted from original code provided by Dr. David Fried %%%%

%%%% Error checking for correct input sizes and enough inputs

%%%Inputs
% dpx - Hartmann sensor phase difference for x-direction

% dpy - Hartman sensor phase difference for y-direction

% sigsqx - Noise variances associated with dpx-values

% sigsy - Noise variances assoiciated with dpy-values

% mask - Array of 1's and 0's, where the 1's in the matrix represent the
% positions that are within the aperture of the sensor

%%%Outputs

% E - Complex phasor representing the reconstructed wave function

% phi - Quasi smooth phase corresponding to E, formed to deal with branch
% cuts in near minimum intensity locations on the grid

% VLQ - Variance like quantity characterizing the variability of E and of
% phi


[a,b]=size(dpx);
if a~=b; error('dpx is not a square array.');end

N=round(log2(a));
if a~=2^N; error('Size of dpx is not 2^N-by-2^N.');end

[c,d]=size(dpy);
if c~=a | d~=a; error('Size of dpy is different from that of dpx.');end

[c,d]=size(sigsqx);
if c~=a | d~=a; error('Size of sigsqx is different from that of dpx.');end

[c,d]=size(sigsqy);
if c~=a | d~=a; error('Size of sigsqy is different from that of dpy.');end

if nargin<4
    error('Not enough inputs. You atleast need X&Y phasers and X&Y noise'); 
end

if nargin==4
    mask=ones(a+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dpu=zeros(a+1,a);
dpv=zeros(a,a+1);
sigsqu=inf*ones(a+1,a);
sigsqv=inf*ones(a,a+1);

%assigning values from Shack Hartmann Sensor to first Hudgin lattice for exponential reconstructor.

for x=1:a
    for y=1:a
        if x+y==2*round((x+y)/2)
            u=(x-y+a)/2+1;
            v=(x+y)/2;
            dpv(v,u)=dpx(y,x)+dpy(y,x);
            sigsqv(v,u)=sigsqx(y,x)+sigsqy(y,x);
        else
            u=(x-y+a+1)/2;
            v=(x+y+1)/2;
            dpu(v,u)=dpx(y,x)-dpy(y,x);
            sigsqu(v,u)=sigsqx(y,x)+sigsqy(y,x);
        end
    end
end

%Hudgin lattice generated for phase 1. Passed into NWRCER with differential phasors calculated in the equation (exp(i*dp#))

[Ep,sigsqp] = NVWCER4SID(exp(i*dpu),exp(i*dpv),sigsqu,sigsqv); %Fitting of data from SHWFS in a Fried/Southwell geometry to Hudgin Geometry

%Create raw lattices so that only values which correpsond to the orignal
%lattice are taken from the hudgin geometry generation and calculation

E1=NaN*ones(a+1,a+1);     
sigsq1=NaN*ones(a+1,a+1);

%Go through point-by-point on large (hudgin/phase 1 lattice) and assign 
%only values which correspond to values overlapping with the original
%lattice. If they are outside the appearture they maintain the value of 0 phase and Inf variance (noise)


for u=1:a+1     
    for v=1:a+1
        x=u+v-1-a/2;
        y=-u+v+1+a/2;
        if ((x>=1 && x<=a+1) & (y>=1 & y<=a+1)) & isfinite(sigsqp(u,v))
            E1(y,x)=Ep(v,u);
            sigsq1(y,x)=sigsqp(v,u);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Assign Hartman data to second lattice and carry out reconstruction%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second lattice reconstruction algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpu=zeros(a+1,a);
dpv=zeros(a,a+1);
sigsqu=inf*ones(a+1,a);
sigsqv=inf*ones(a,a+1);

for x=1:a
    for y=1:a
        if x+y==2*round((x+y)/2)
            u=(x-y+a)/2;
            v=(x+y)/2;
            dpu(v,u)=dpx(y,x)-dpy(y,x);
            sigsqu(v,u)=sigsqx(y,x)+sigsqy(y,x);
else
            u=(x-y+a+1)/2;
            v=(x+y-1)/2;
            dpv(v,u)=dpx(y,x)+dpy(y,x);
            sigsqv(v,u)=sigsqx(y,x)+sigsqy(y,x);
        end
    end
end

[Epp,sigsqpp] = NVWCER4SID(exp(i*dpu),exp(i*dpv),sigsqu,sigsqv);

E2=NaN*ones(a+1,a+1);
sigsq2=NaN*ones(a+1,a+1);

for u=1:a+1
    for v=1:a+1
        x=u+v-a/2;
        y=-u+v+1+a/2;
        if ((x>=1 & x<=a+1) & (y>=1 & y<=a+1)) & isfinite(sigsqp(u,v)) 
            E2(y,x)=Epp(v,u);
            sigsq2(y,x)=sigsqpp(v,u);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MERGE THE LATTICE MATRICES BACK TOGETHER INTO ORIGINAL%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E2t=NaN*ones(a+1);
sigsq2t=NaN*ones(a+1);

for x=1:a+1
    for y=1:a+1
        if x+y~=2*round((x+y)/2)
            S=0; T=0;
            if x-1>=1
                S=E1(y,x-1)/sigsq1(y,x-1);
                T=1/sigsq1(y,x-1);
            end
            if x+1<=a+1
                S=S+E1(y,x+1)/sigsq1(y,x+1);
                T=T+1/sigsq1(y,x+1);
            end
            if y-1>=1
                S=S+E1(y-1,x)/sigsq1(y-1,x);
                T=T+1/sigsq1(y-1,x);
            end
            if y+1<=a+1
                S=S+E1(y+1,x)/sigsq1(y+1,x);
                T=T+1/sigsq1(y+1,x);
            end
            T=4/T;
            S=S/(4/T);
            E2t(y,x)=S;
            sigsq2t(y,x)=T;
        end
    end
end
    

% As an option the estimated amplitude of the reconstructed wave funcion at
% each of the corners of the Hartmann sensor's sub apertures can be provided
% as an input. If this input is provided the amplitudes of the reconstructed
% wave function will be adjusted to match these values. If not provided then
% the reconstructed wave function will have unity amplitude at all points.
% 
% The input called mask indicates the positions corresponding to points that
% are within the aperture. The elements (positions) of mask correspond to
% corners of the Hartmann sensor's sub apertures, for each of the elements of
% which one of the elements of dpx, dpy, sigsqx,
% and sigsqy correspond. Any element of sigsqx and of sigsqy which has a
% value of infinity corresponds to a sub aperture for which the four corners
% all have mask element values of zero.

[x,y]=meshgrid(1:a+1);
ii=find(x+y~=2*round((x+y)/2)); 
S=sum(E2t(ii).*conj(E2(ii)).*mask(ii)./(sigsq2t(ii)+sigsq2(ii))); 
phi21=angle(S);
E=E1;
E(ii)=E2(ii)*exp(i*phi21);
VLQ=sigsq1;
VLQ(ii)=sigsq2(ii);
jj=find(~isfinite(VLQ));
E(jj)=NaN;

E=E./abs(E);
if nargin==6
    E=amp.*E;
end
phi=BranchPointPhase(E);
end


% BranchPointPhase generates a phase function, PHI that corresponds to
% |2*pi| to the phase of the input complex phasor, U, and that is as far as
% possible smooth. i.e, minimisez differences between values of phi between
% adjacent points in the array of values of phi. The magnitude of the
% difference of phase values between adjacent points should be less than pi
% everywhere except branch-cuts.

% The input function, U, defines the extent of the aperture by its value of
% being equal to NaN for those array points corresponding to positions
% outside the aperture. The U input array is a square array and is smaller
% than the coresponding circular aperture. 

% PHI is calculated by evaluating the difference of the phase between
% adjecent values of U, using phase differences dPx and dPy, to calculate
% the curl for each elemental square of the array, and using the values of
% the curl to determine the lcoation and sign of the branch points. Using
% this information about the branch points, the complex phasor Uh,
% associated with the branch points is calculated and from this the hidden
% phase is phi, is calculated as a result. The ratio U over Uh is denoted
% by n and is calculated as well. Then the phase differences dpx and dpy
% associated with u are evaluated. 
% The associated phase is generated by adding/subtracting dpy values
% alongside the central line of the array to get the values of all the
% points above/below the center point of the central line. Then the dpx
% values from the central line for each point to the right/left of the
% central line. The output phase values, PHI is a result of the adding the
% hidden phases with the array resulted by the calculation of associated
% phase.

% If the number of possitive sign branch points is not equal to the
% negative ones (missing, or outside the aperture), an additional branch
% point is provided. This point is placed outside the aperture and at a
% location chosen on the basis of a potential field, V, formed by the
% possitive and negative sign branch points. 

% Inputs: U= complex phasor array (if outside the aperture, NaN values in
% the square matrix)
%
%
% Outputs: PHI (N,M) = Maximally smooth phase corresponding to modulo 2*pi
% to the phase of U (values outside the aperture are NaN).
%
%          BPcount (1,1) = Number of branch points in U
%          phi (N,M) = Hidden phase corresponding to U
%
%          BPes (a, 3) = Second column of values gives the x-axis,
%          possition of branch point (column number) while the first column 
%          gives the y axis (row number). And thir column indicates by +/-1
%          if the branch point is negative or possitive

function [PHI,BPcount,phi,BPes]=BranchPointPhase(U)
[N,M]=size(U);
if N~=M; error('Array should be square.'); end 

%Equations 35a/b
%Atan2 is for arg() which calculates principle value form

dPx=U(:,2:end).*conj(U(:,1:end-1));
dPx=atan2(imag(dPx),real(dPx));
dPy=U(2:end,:).*conj(U(1:end-1,:));
dPy=atan2(imag(dPy),real(dPy));

BPes=[]; 
curl=dPx(1:end-1,:)+dPy(:,2:end)-dPx(2:end,:)-dPy(:,1:end-1); %Eqn 36 for n=1:N-1
for n=1:N-1
    for m=1:M-1
        if abs(curl(n,m))>.1*pi
            BPes=[BPes; n m sign(curl(n,m))];
        end
    end
end

BPcount=size(BPes,1);
bpcount=BPcount;
Uh=ones(N,M);

if bpcount>0
    
    BPexcess=sum(BPes(:,3));
    %Concerned with ensuring the number of positive and negative branch
    %points are the same. And if not carry out this process
    
    while BPexcess~=0               %Calculating Step 2 1/2 (eqn 44)
        R=(N+3)/2;
        theta=(0:359)*pi/180;
        x=R*cos(theta)+M/2;
        y=R*sin(theta)+N/2;
        V=zeros(1,360);
        for k=1:BPcount
            V=V+BPes(k,3)./sqrt((x-BPes(k,2)).^2+(y-BPes(k,1)).^2);
        end
        [mx,ii]=max(BPexcess*V);
        bpcount=bpcount+1;
        BPes=[BPes; y(ii) x(ii) -sign(BPexcess)];
        BPexcess=sum(BPes(:,3));
    end
    [x,y]=meshgrid(1:M,1:N);
    for bpc=1:bpcount
        X=BPes(bpc,2)+0.5;
        Y=BPes(bpc,1)+0.5;
        pn=BPes(bpc,3);
        if pn>0
            Uh=Uh.*[(x-X)+1i*(y-Y)];
        else
            Uh=Uh./[(x-X)+1i*(y-Y)];
        end
    end
end

phi=angle(Uh);
ii=find(~isfinite(U));
phi(ii)=NaN;
u=U./Uh;
dpx=u(:,2:end).*conj(u(:,1:end-1));
dpx=atan2(imag(dpx),real(dpx));
dpy=u(2:end,:).*conj(u(1:end-1,:));
dpy=atan2(imag(dpy),real(dpy));
Phi=Reconstructor(dpx,dpy);
PHI=Phi+phi;
end

%
% Function Reconstructor, used to accomplish reconstruction based on simple
% addition of phase differences. Starting from center of array, first the
% phase differences along y-axis are added to central line. Then starting
% from each point in that line the phase differences are added along the
% x-asi direction. 
%
% Central element is always 0, the mean phase of all the elements will be
% 0.
function phi=Reconstructor(pdx,pdy)

if nargin==2
    yn=-1;
end

[a,b]=size(pdx);
if a~=b+1
    error('Size of pdx is not N-by-(N-1).')
end
[c,d]=size(pdy);
if (d~=a) | (c~=b)
    error('Size of pdy does not properly correspond to that of pdx.')
end

N=size(pdx,1);
hN=round(N/2);
pd=-flipud(pdy(1:hN-1,hN));
phi=flipud(cumsum(pd));
pd=pdy(hN:N-1,hN);
phi=[phi; 0; cumsum(pd)];
pd=-fliplr(pdx(:,1:hN-1));
phi=fliplr(cumsum([phi pd],2));
pd=pdx(:,hN:N-1);
phi=[phi(:,1:hN-1) cumsum([phi(:,hN) pd],2)]; %flip array left to right (added ;)
W=isfinite(phi);
phi=phi-mean(phi(W(:)));
phi(~W)=0;

end