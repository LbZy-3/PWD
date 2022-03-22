 clear all
 close all


%  file=uigetfile('*.png');
%  B=imread(file);
%  sig=rgb2gray(B);

u=1;             % length unit, e.g. 1 µm




lambda=0.5*u;

x=[-5:0.005:5]*u
y=[-5:0.005:5]*u
Nx=size(x,2);
Ny=size(y,2);

A=zeros(Ny,Nx);
A(floor(Ny/2),floor(Nx/2))=1;
sig=A;

freq_array1=(0:Nx-1)-Nx/2; %possible spatial freq. in the signal (no need of x)
freq_array2=freq_array1/(max(x)-min(x)); % spatial freq. in x unit
kx=2*pi*freq_array2;  % angular spatial freq in rad/x_unit  (kx)
freq_array1=(0:Ny-1)-Ny/2; %possible spatial freq. in the signal (no need of y)
freq_array2=freq_array1/(max(y)-min(y)); % spatial freq. in y unit
ky=2*pi*freq_array2;  % angular spatial freq in rad/x_unit  (ky)

w = 1;%window2(Ny, Nx, @gausswin); 

TF=fftshift(ifft2((sig))).*w; % centred or.
SIG=(fft2((TF))); %idem
%-------------------------------------------------------
% Definit° for Transfer function H=exp(ikz*z) & propagation
%---------------------------------------------------------
k=2*pi/lambda;
[KX, KY] = meshgrid(kx,ky);
KZ=(k.^2-KX.^2-KY.^2).^0.5; % Array of kz


Hat=(imag(KZ)==0);

%-------------------------------------------------------
% DISPLAY SIGNAL AND FFT in object plane
%-------------------------------------------------------
% SET "if 1" to activate display, "if 0" to deactivate

if 1 % 
    figure(1)
    subplot(3,1,1)
    imagesc(x,y,sig);
    title('sig');
    subplot(3,1,2)
    imagesc(kx,ky,abs(TF+Hat));
    title('abs(fft(sig))');
    subplot(3,1,3)
    imagesc(kx,ky,angle(TF));
    title('angle(fft(sig))');

end


zz=[0.02:0.1:1]*u; %range of propagation with step increment
z=0;
figure(2);
for step=[1:1:size(zz,2)]; % forward propagation steps (arbitrary) in x unit
    z=zz(step) % forward propagation steps (arbitrary) in x unit      
    H=exp(1i*KZ*(zz(step))); % Transfer function
    TFz=TF.*H;
    im=(fft2(TFz)); 
    imagesc(x,y,abs(im));
    colormap gray
    %TF=fftshift(ifft2((im)));   
    %z=zz(step);
    title(['Propagation to z=' num2str(zz(step))]);  
    pause(0.2)
end
TF=TF.*Hat;
for step=[size(zz,2):-1:1]; % forward propagation steps (arbitrary) in x unit
    z=zz(step) % backward propagation steps (arbitrary) in x unit      
    H=exp(1i*KZ*(zz(step))); % Transfer function
    TFz=TF.*H;
    im=(fft2(TFz)); 
    imagesc(x,y,abs(im));
    colormap gray
    %TF=fftshift(ifft2((im)));   
    %z=zz(step);
    title(['Propagation to z=' num2str(zz(step))]);  
    pause(0.2)
end
