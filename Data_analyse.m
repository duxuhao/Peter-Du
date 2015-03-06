%%for analyse the data for single 
%%read data and set parameters
WhiteNoise=wavread('noise.wav');
SM=wavread('SM_2_1.wav');
fs=44100;%recorded sample frequency
fs2=10000;%the orginal signal frequency
p0=2e-5;%reference pressure
%%
%calculate the spectrum of the origin signal
fftWhiteNoise=abs(fft(WhiteNoise,fs2)*2)/fs2;
ReSM=resample(SM,fs2,fs);%resample to fit the original signal
SMNoise=ReSM((end-fs2+1):end);
fftSMNoise=abs(fft(SMNoise,fs2)*2)/fs2;
figure
plot(20*log10(fftSMNoise),'r');%test for the edit
xlabel('Frequency (Hz)','fontsize',14);
ylabel('dB','fontsize',14);
axis([0 fs2/2 min(20*log10(fftSMNoise)) max(20*log10(fftSMNoise))])
%%
%use the comb-filter to process the data
f0=50;%
beta=-0.3;
data=ReSM;
FS=fs2;%the frequency we use in the calcualtion
%SM2=My_Comb(ReSM,fs2,f0,beta);%SM2 the signal after the comb filter
SM2=My_Comb(data,FS,f0,beta);
%SM2=ReSM;
%FS=fs;
SM2Noise=SM2((end-FS+1):end);
L=length(SM2);
SM2Signal=data((round(L/2)-FS+1):round(L/2));
fftSM2Noise=abs(fft(SM2Noise,FS)*2)/FS;
fftSM2Signal=abs(fft(SM2Signal,FS)*2)/FS;
hold on
plot(20*log10(fftSM2Noise));
legend('Origin Noise','Measured Noise after Comb filter',0)
title('Spectrum of the Measured Noise','fontsize',14)
set(gca,'fontsize',14)
%%
%calculate the transfer function
fftSM2difference=fftSM2Signal.^2-fftSM2Noise.^2;
for n=2:FS;
    if fftSM2difference(n)<=0
        fftSM2difference(n)=fftSM2difference(n-1);
    end
end
fftSM2difference2=sqrt(fftSM2difference);
FrequencyResponse=fftSM2difference2./fftWhiteNoise;
ImpulseResponse=ifft(FrequencyResponse,fs2);
%draw the figure
%{
figure
subplot(2,2,1)
plot(20*log10(abs(fftSM2difference/p0)),'k','linewidth',2);
xlabel('Frequency (Hz)');
axis([0 fs2/2 min(20*log10(abs(fftSM2difference/p0)))-10 max(20*log10(abs(fftSM2difference/p0)))+10])
title('Difference between Measured Signal and Noise')
subplot(2,2,2)
plot(20*log10(abs(FrequencyResponse/p0)),'k','linewidth',2);
xlabel('Frequency (Hz)');
axis([0 fs2/2 min(20*log10(abs(FrequencyResponse/p0)))-10 max(20*log10(abs(FrequencyResponse/p0)))+10])
title('Transfer Function')
subplot(2,2,3)
plot(abs(ImpulseResponse),'k');
xlabel('Frequency (Hz)');
axis([0 fs2/2 0 max(abs(ImpulseResponse))])
title('Impulse Response')
subplot(2,2,4)
plot((1:length(ReSM))/fs2,20*log10(abs(ReSM)),'r')
hold on;
plot((1:length(ReSM))/fs2,20*log10(abs(SM2)));
xlabel('Time (s)');
title('Time Sequence')
legend('Origin Signal','Signal after Comb filter',0)
%}
%%
%FS=fs;
figure
subplot(2,2,1)
%calculate the belonged Band Number
ExtractFrequency=2000;
Start=round(8.9*FS);
End=round(9.2*FS);
OctaveBand=[31.5 63 125 250 500 1000 2000 4000];
OneThirdOctaveBand=[20 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000];
for i=1:length(OneThirdOctaveBand)%find the coresponse band
    if ExtractFrequency < OneThirdOctaveBand(i)*(2^(1/6)) && ExtractFrequency > OneThirdOctaveBand(i)/(2^(1/6))
        ExtractFrequency=OneThirdOctaveBand(i);break
    end
end
bandwidth=round(ExtractFrequency*(2^(1/6))-ExtractFrequency/(2^(1/6)));%set to the bandpass filter's bandwidth
%bandwidth=10;
D = round(bandwidth/i);%the edge decay
%D=5;
%SM_Filte=ExtractFre(SM,ExtractFrequency,bandwidth,FS);%use the origin signal to
%calculate the decay time
SM_Filte=ExtractFre(SM2,ExtractFrequency,bandwidth/2,FS,D);
title(['Designed Filter Curve ',num2str(ExtractFrequency),'Hz and ' ,num2str(bandwidth),' Hz Bandwidth'],'fontsize',14);
%legend([num2str(bandwidth),' Hz Bandwidth'],'Location','South')
subplot(2,2,2)
plot(20*log10(abs(fft(SM_Filte,FS)*2/FS)),'k','linewidth',2)
legend(['Edge Decay ',num2str(D),' dB'],0)
title(['Spectrum(',num2str(ExtractFrequency),'Hz Filter)'],'fontsize',14)
ylabel('dB','fontsize',14);
xlabel('Frequency (Hz)','fontsize',14);set(gca,'fontsize',14)
axis([0 FS/2 min(20*log10(abs(fft(SM_Filte,FS)*2/FS)))-10 max(20*log10(abs(fft(SM_Filte,FS)*2/FS)))+10])
%figure
subplot(2,2,3)
plot((1:length(SM_Filte))/FS,20*log10(abs(SM_Filte)),'k')
ylabel('dB','fontsize',14);
xlabel('Time (s)','fontsize',14);set(gca,'fontsize',14)
axis([Start/FS-0.1 End/FS+0.1 -150 0]);
title(['Time Sequence (',num2str(ExtractFrequency),'Hz Filter)',],'fontsize',14)
%Start=input('start Number:');
%End=input('End Number:');
subplot(2,2,4)
figure
yy=My_BackwardIntergration(SM_Filte,Start,End,FS);
plot((Start:End-1)/FS,10*log10(yy)-10*log10(yy(1)),'k','linewidth',2);
ylabel('dB','fontsize',14);grid on;
xlabel('Time (s)','fontsize',14);set(gca,'fontsize',14)
axis([Start/FS (End+2)/FS min(10*log10(yy)-10*log10(yy(1))) 0])
%axis([1/FS (End-Start+2)/FS min(10*log10(yy)-10*log10(yy(1))) 0])
title([num2str(ExtractFrequency),' Hz Decay Sequence'],'fontsize',14)