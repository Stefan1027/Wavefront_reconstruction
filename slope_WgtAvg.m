function [Fx,Fy,Magnitudes,SNR] = slope_WgtAvg(phase2pi,amplitude2pi,NN,SNR0,flag)

%Create Complex Matrix from amplitude and phase

Phi_Complex=amplitude2pi.*exp(1i*phase2pi);

%Generate complex differentials

PhiDx=Phi_Complex(:,2:end).*conj(Phi_Complex(:,1:end-1));
PhiDy=Phi_Complex(2:end,:).*conj(Phi_Complex(1:end-1,:));

%Magnitude of differentials for weighting

MagX=abs(PhiDx);
MagY=abs(PhiDy);

%Phase of difference

thetaX=angle(PhiDx);
thetaY=angle(PhiDy);
Intensity=abs(Phi_Complex).^2;
Intensity=[Intensity(:,1) Intensity Intensity(:,end)];
Intensity=[Intensity(1,:); Intensity; Intensity(end,:)];

%Repeat edge values

MagX=[MagX(:,1) MagX MagX(:,end)];
MagY=[MagY(1,:); MagY; MagY(end,:)];
thetaX=[thetaX(:,1) thetaX thetaX(:,end)];
thetaY=[thetaY(1,:); thetaY; thetaY(end,:)];

%Calculate differential size

[y,x]=size(phase2pi);
dy=round(y/NN);
dx=round(x/NN);

%Slope outputs

Fx=zeros(NN);
Fy=zeros(NN);
Magnitudes=zeros(NN);

%Weighting Factors for matrix manipulation since the column (row) for x(y) 
%differentials has overlapping of sample area for each lenslet and is 
%ideally located in the middle of each phase point (fried geometry) a 
%weight of one half is applied to the over sampled values on the edges

Wx=ones(dy,dx+1);
Wx(:,1)=0.5;
Wx(:,end)=0.5;
Wy=ones(dy+1,dx);
Wy(1,:)=0.5;
Wy(end,:)=0.5;
Wp=ones(dy+2,dx+1);
Wp(1,:) = 0.5;
Wp(end,:)=0.5;
WF=ones(dx+2, dx+2)
WF(1,:)=0.2;
WF(end,:)=0.2;
Wmag=ones(dy+1,dx+1);
Wmag(:,1)=0.5;
Wmag(:,end)=0.5;
Wmag(1,:)=Wmag(1,:).*.5;
Wmag(end,:)=Wmag(end,:).*.5;
thetaFx=zeros(dy,dx+1);
thetaFy=zeros(dy+1,dx);
magFx=zeros(dy,dx+1);
magFy=zeros(dy+1,dx);

for u=1:NN
    for v=1:NN
        thetaFx=thetaX((v-1)*dy+1:v*dy,(u-1)*dx+1:u*dx+1);
        thetaFy=thetaY((v-1)*dy+1:v*dy+1,(u-1)*dx+1:u*dx);
        magFx=MagX((v-1)*dy+1:v*dy,(u-1)*dx+1:u*dx+1).*Wx;
        magFy=MagY((v-1)*dy+1:v*dy+1,(u-1)*dx+1:u*dx).*Wy;
        Fx(v,u)=sum(sum(thetaFx.*magFx))/sum(sum(magFx)).*dx; %tiltx
        Fy(v,u)=sum(sum(thetaFy.*magFy))/sum(sum(magFy)).*dy; %tilty
        Magnitudes(v,u)=sum(sum(Intensity((v-1)*dy+1:v*dy+1,(u-1)*dx+1:u*dx+1).*Wmag));
    end
end

SNR=SNR0*Magnitudes/mean(Magnitudes(:));
delX=randn(NN,NN).*pi./SNR;
delY=randn(NN,NN).*pi./SNR;

if flag==1
    Fx=Fx+delX;
    Fy=Fy+delY;
end

end
