function [Fx,Fy,Magnitudes,SNR,NN] = wavefront_sensor(deviation_x,deviation_y,NN_sensor,SNR0,flag,spot_intensity)
f=4.1; %effective focal length, however needs to be outputed more accurately from the WFS, to increase accuracy
Fx=deviation_x/f;
Fy=deviation_y/f;
NN=NN_sensor;
Magnitudes=spot_intensity;

SNR=SNR0*Magnitudes/mean(Magnitudes(:));
delX=randn(NN,NN).*pi./SNR;
delY=randn(NN,NN).*pi./SNR;

if flag==1
    Fx=Fx+delX;
    Fy=Fy+delY;
end

end