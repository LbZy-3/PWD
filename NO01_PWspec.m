% Simulate Forward or backward propagation from near field signal
close all,clear all,clc
%load near-field signal:
[x,sig]=nearfield_signal(1);  % sig: plotted later on
N=length(x); 
L=1;% wavelength in unit of 'x', e.g. 1µm

% Array of frequencies centered on 0, several options:
freq_array1=(0:N-1)-N/2; %possible spatial freq. in the signal (no need of x)
freq_array2=freq_array1/(max(x)-min(x)); % spatial freq. in x unit
kx=2*pi*freq_array2;  % angular spatial freq in rad/x_unit  (kx)
f=kx;  %<= frquency choice in absica


% fft computation
%% 1)Simplest expression:  the spectrum have 0th freq at first point
% TF=fft(sig);        %can't be used here later on, cause.. 
% SIG=ifft((TF));     %.. we have centered our freq in middle
%% (2) More common expression: the spectrum have 0th freq in middle
%TF=fftshift(fft(sig));
%SIG=ifft(ifftshift(TF));
%% (3) like 2) but phase Origin @center, (not @left as both 2 above)
TF=fftshift(fft(ifftshift(sig))); % centred or.
SIG=fftshift(ifft(ifftshift(TF))); %idem

%-------------------------------------------------------
% DISPLAY SIGNAL AND FFT in object plane
%-------------------------------------------------------
% SET "if 1" to activate display, "if 0" to deactivate
if 1 % 
subplot(3,1,1)
 plot(x,sig,'--ro','markersize',1,'linewidth',2),ylabel('sig'); axis tight;
 title(['Signal with following number of points: ' num2str(N)])
 subplot(3,1,2)
 plot(f,abs(TF)); axis tight; ylabel('abs(fft(sig))')
 subplot(3,1,3)
 plot(f,angle(TF),':ro','markersize',1);  axis tight;
 ylabel('angle(fft(sig))'), legend('Phase')  
 xlabel('spatial freq. (if angular spatial freq: rad/µm)')
pause();clf
end 

 %-------------------------------------------------------
 % Definit° for Transfer function H=exp(ikz*z) & propagation
 %---------------------------------------------------------
k=2*pi/L;
kz=(k.^2-kx.^2).^0.5; % Array of kz

zz=[0.02:0.1:1]; %range of propagationwith step increment

for step=[1:1:10]; % forward propagation steps (arbitrary) in x unit
z=zz(step)
H=exp(1i*kz*z); % Transfer function
TFz=TF.*H;
SIGz(step,:)=fftshift(ifft(ifftshift(TFz))); % <=use if choice (3)
%SIGz=ifft(ifftshift(TFz));          % <=use if choice (2) 
plot(x,abs(SIGz));xlabel('x');title(['Propagation to z=' num2str(z)])
pause()
end
figure
imagesc(x,-zz,log(abs(SIGz)+0.001));axis off; 
title(['From 0 to z=' num2str(z) ' Log scale'])
colormap hot
pause
disp('Removing evanescent waves & refocus')
pause
%--------------------
% REmoving evanescent components from TF
%-----------------
plot(kx,real(kz)),xlabel('k_x');title('real(kz)');pause
plot(kx,imag(kz)),xlabel('k_x');title('imag(kz)');pause
Hat=(imag(kz)==0); % plot(kx,Hat); %to check
plot(kx,Hat);xlabel('k_x');title('Filtre (hat)');pause
TF=TF.*Hat;


% Backpropagation
clear SIGz
for step=[length(zz):-1:1]; 
z=zz(step)
H=exp(1i*kz*z); % Transfer function
TFz=TF.*H;
SIGz(step,:)=fftshift(ifft(ifftshift(TFz))); % <=use if choice (3)
%SIGz=ifft(ifftshift(TFz));          % <=use if choice (2) 
plot(x,abs(SIGz));xlabel('x');title(['Propagation to z=' num2str(z)])
pause
end
figure
imagesc(x,-zz,log(abs(SIGz)+0.001));axis on; 
title(['Refocusing.. (Log scale)'])
colormap hot


 disp('The end')