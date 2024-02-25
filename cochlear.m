clear all;close all; clc;
fs=48000;
f=2000;
n=0.01:0.01:100;
y=sin(2*pi*(f/fs)*n );
[y, fs]=audioread('Music.wav');
figure;spectrogram(y);
T=1/fs;
%% Bandpass Filter Design Parameters for mid and outer ear
% Bandpass Filter Design Parameters for outer ear

cf = 3000; 
bw = 4000; 
fo = 4;
dg = 20; 
dgl = 10^(dg / 20);

ncf = cf / (fs / 2);
nbw = bw / (fs / 2);

% co-efficient of filter using butterworth filter
[bo, ao] = butter(fo, [ncf - nbw/2 ncf + nbw/2], 'bandpass');
bo = bo * dgl;

figure; freqz(bo, ao,  fs); xlabel('Frequency (Hz)');ylabel('Magnitude');title('Bandpass Filter Frequency Response(OE)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandpass Filter Design Parameters for middle ear
cf = 3000; 
bw = 1000; 
fo = 6;
dg = 20; 
dgl = 10^(dg / 20);

ncf = cf / (fs / 2);
nbw = bw / (fs / 2);

% co-efficient of filter using butterworth filter
[bm, am] = butter(fo, [ncf - nbw/2 ncf + nbw/2], 'bandpass');
bm = bm * dgl;

figure;freqz(bm, am,  fs); xlabel('Frequency (Hz)');ylabel('Magnitude');title('Bandpass Filter Frequency Response(ME)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cascaded filter of outer and middle ear
b=conv(bo, bm);
a=conv(ao, am);
y=filter(b,a,y);
figure; freqz(b,a,fs);title('outer+middle ear filter');
figure; zplane(b,a);
sys=filt(b,a);

%% inner ear (Basilar Membrance output)
div =128;d=3.5/div;
qp = linspace(10, 5.5, div);
qz = linspace(22, 12, div);

fp=zeros(1,div);
fz=zeros(1,div);
bwp=zeros(1, 128);
bwz=zeros(1, 128);
for i=1:div
    fp(i)=20000*10^(-(0.667*(i-1)*d));
    fz(i)=fp(i)*1.0429;
    bwp(i)=fp(i)/qp(i);
    bwz(i)=fz(i)/qz(i);
end
fil_coff=zeros(128, 9);
% k=fil_coff(i,1) 
% g0=fil_coff(i,2) = 1-a0
% a0=fil_coff(i,3) 
% gp=fil_coff(i,4) = 1-b1+b2
% b1=fil_coff(i,5) 
% b2=fil_coff(i,6) 
% gz=fil_coff(i,7) = 1/(1-a1+a2)
% a1=fil_coff(i,8) 
% a2=fil_coff(i,9)
for i=1:div
    p1=bwp(i)/(2*qz(i));q1=p1*sqrt(4*qz(i)^2-1);
    fil_coff(i,8)=2*exp(-p1*T ).*cos(q1*T );fil_coff(i,9)=exp(-2*p1*T );

    p2=bwp(i)/(2*qp(i));q2=p2*sqrt(4*qp(i)^2-1);
    fil_coff(i,5)=2*exp(-p2*T )*cos(q2*T );fil_coff(i,6)=exp(-2*p2*T );
    
    flp=fz(i)*1.4;
    thetalp=2*pi*(flp/fs);

    fil_coff(i,3)=2-cos(thetalp)-sqrt((2-cos(thetalp))^2-1);
    
    fil_coff(i,1)=fp(i)/fz(i);
a0=fil_coff(i,3) ;
    fil_coff(i,2)=1-a0;
b1=fil_coff(i,5) ;
b2=fil_coff(i,6) ;
    fil_coff(i,4)=1-b1+b2;
a1=fil_coff(i,8) ;
a2=fil_coff(i,9) ;
    fil_coff(i,7)=(1/1-a1+a2);
num=[1-a0];
den=[1 -a0];
    y=filter(num,den,y);

num=[1-b1-b2];
den=[1 -b1 -b2];
    y=filter(num,den,y);
    dis_o(:,i)=y;

num=[1 -a1 -a2];
den=[1-a1-a2];
    y=filter(num,den,y);
    pre_o(:,i)=y;
end

%% differentiation
d_dis_o=zeros(size(dis_o)); d_pre_o=zeros(size(pre_o));
d_dis_o(:, 1)=dis_o(:,1);
d_pre_o(:, 1)=pre_o(:,1);
for i=2:div
    d_dis_o(:, i)=dis_o(:,i)-dis_o(:, i-1);
    d_pre_o(:, i)=pre_o(:,i)-pre_o(:, i-1);
end
for i=2:length(dis_o)
    d_dis_o(i,:)=dis_o(i, :)-dis_o(i-1,:);
    d_pre_o(i,:)=pre_o(i, :)-pre_o(i-1,:);
end

clear a0;clear num;clear den; clear flp; clear thetalp;
clear b1;clear b2;clear a1;clear a2;
clear p1;clear p2;clear q1;clear q2;

h_o=zeros(size(dis_o));
%% hair cell
for i=1:128
    c0=exp(-30*(2*pi/48000));
    num=[1-c0];
    den=[1 -c0];
    h_o(:,i)=filter(num, den, dis_o(:,i));
end

plot(h_o(:,81))