clear all;

%Sampling points
Np = 1e6;

%Number of the grating layers
N = 50;
c = 3e8;
lamda0 = 14.1e-6;
L = N * lamda0;

%Set the modal gain g0
Lc = 5e-3;
g0 = 1 / Lc ;
alpha_h = 4;

%Frequency initialisation
%Set the range of frequency
fs = 0; fe = 2; 
%Set the FWHM - full width at half maximum
FWHM_f = 300e9;
%Centre frequency in THz
fB=2.89;

frequency = fs*fB:(fe-fs)*fB/Np:fe*fB;

%Set the refractive indices
n1=3.717; n2=3.5;
n_offset = 3.605;

%In the case of non-grating strucuture 
n1_0=3.717; n2_0=n1_0;

%Define tao_0
tao_0 = N / (2 * fB);

%Set the gain coefficient g (g_frequency = 1/2 g_power)
gL = 0;
g = gL / 200;

%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];

S12_0=(1/(2*n1_0))*[n1_0+n2_0 n1_0-n2_0;n1_0-n2_0 n1_0+n2_0];
S21_0=(1/(2*n2_0))*[n2_0+n1_0 n2_0-n1_0;n2_0-n1_0 n2_0+n1_0];

for i = 1:length(frequency)                                              
fn=frequency(i)/fB;    
%Calculate the optical thickness
theta1 = (pi/2)*fn+1i*g;       %For slot     
theta2 = (pi/2)*fn+1i*g;       %For lamda      

theta2_1=(pi/2)*fn+1i*g;       %For lamda         
theta2_2=(pi)*fn+1i*2*g;       %For a quarter-wave phase-shift

theta3_1 = (pi/2)*fn+1i*g;          %For lamda       
theta3_2 = (3*pi/2)*fn+1i*3*g;      %For 2*lamda     
theta3_3 = (3*pi)*fn+1i*6*g;        %For 3.5*lamda           
theta3_4 = (2*pi)*fn+1i*4*g;        %For 2.5*lamda 
%Calculate the medium constant P
Pn1 = [exp(1i*theta1) 0;0 exp(-1i*theta1)];
Pn2 = [exp(1i*theta2) 0;0 exp(-1i*theta2)];

Pn2_1=[exp(1i*theta2_1) 0;0 exp(-1i*theta2_1)];
Pn2_2=[exp(1i*theta2_2) 0;0 exp(-1i*theta2_2)];

Pn3_1=[exp(1i*theta3_1) 0;0 exp(-1i*theta3_1)];
Pn3_2=[exp(1i*theta3_2) 0;0 exp(-1i*theta3_2)];
Pn3_3=[exp(1i*theta3_3) 0;0 exp(-1i*theta3_3)];
Pn3_4=[exp(1i*theta3_4) 0;0 exp(-1i*theta3_4)];
%M matrix for the periodic grating strcture
M1 = (Pn1*S12*Pn2*S21)^N;

M2 = ((Pn1*S12*Pn2_1*S21)^24) * (Pn1*S12*Pn2_2*S21) * ((Pn1*S12*Pn2_1*S21)^25);

M3 = ((Pn1*S12*Pn3_1*S21)^13) * (Pn1*S12*Pn3_2*S21)*(Pn1*S12*Pn3_1*S21)*(Pn1*S12*Pn3_2*S21)*(Pn1*S12*Pn3_3*S21)*(Pn1*S12*Pn3_2*S21)*(Pn1*S12*Pn3_1*S21)*(Pn1*S12*Pn3_2*S21)*(Pn1*S12*Pn3_1*S21)*(Pn1*S12*Pn3_1*S21)*(Pn1*S12*Pn3_1*S21)*(Pn1*S12*Pn3_4*S21)*(Pn1*S12*Pn3_2*S21)*(Pn1*S12*Pn3_2*S21) * ((Pn1*S12*Pn3_1*S21)^14);
%M matrix for the non-grating structure
M1_0 = (Pn1*S12_0*Pn2*S21_0)^N;

M2_0 = ((Pn1*S12_0*Pn2_1*S21_0)^24) * (Pn1*S12_0*Pn2_2*S21_0) * ((Pn1*S12_0*Pn2_1*S21_0)^25);

M3_0 = ((Pn1*S12_0*Pn3_1*S21_0)^13) * (Pn1*S12_0*Pn3_2*S21_0)*(Pn1*S12_0*Pn3_1*S21_0)*(Pn1*S12_0*Pn3_2*S21_0)*(Pn1*S12_0*Pn3_3*S21_0)*(Pn1*S12_0*Pn3_2*S21_0)*(Pn1*S12_0*Pn3_1*S21_0)*(Pn1*S12_0*Pn3_2*S21_0)*(Pn1*S12_0*Pn3_1*S21_0)*(Pn1*S12_0*Pn3_1*S21_0)*(Pn1*S12_0*Pn3_1*S21_0)*(Pn1*S12_0*Pn3_4*S21_0)*(Pn1*S12_0*Pn3_2*S21_0)*(Pn1*S12_0*Pn3_2*S21_0) * ((Pn1*S12_0*Pn3_1*S21_0)^14);

%Calculate the transmission coefficient
Rtrans1 = 1 / M1(1,1);
Rtrans1_0 = 1 / M1_0(1,1);

Rtrans2 = 1 / M2(1,1);
Rtrans2_0 = 1 / M2_0(1,1);

Rtrans3 = 1 / M3(1,1);
Rtrans3_0 = 1 / M3_0(1,1);
%Calculate the phase 
phase1(i) = angle(Rtrans1);
phase1_0(i) = angle(Rtrans1_0);

phase2(i) = angle(Rtrans2);
phase2_0(i) = angle(Rtrans2_0);

phase3(i) = angle(Rtrans3);
phase3_0(i) = angle(Rtrans3_0);

end

for k = 1:length(frequency)
fn(k)=frequency(k)/fB;   
if k==1
   phase_accu1(k) = abs(phase1(k));
   phase_accu1_0(k) = abs(phase1_0(k));

   phase_accu2(k) = abs(phase2(k));
   phase_accu2_0(k) = abs(phase2_0(k));

   phase_accu3(k) = abs(phase3(k));
   phase_accu3_0(k) = abs(phase3_0(k));
else
   phase_accu1(k) = abs(abs(phase1(k)) - abs(phase1(k-1))) + phase_accu1(k-1);
   phase_accu1_0(k) = abs(abs(phase1_0(k)) - abs(phase1_0(k-1))) + phase_accu1_0(k-1);

   phase_accu2(k) = abs(abs(phase2(k)) - abs(phase2(k-1))) + phase_accu2(k-1);
   phase_accu2_0(k) = abs(abs(phase2_0(k)) - abs(phase2_0(k-1))) + phase_accu2_0(k-1);

   phase_accu3(k) = abs(abs(phase3(k)) - abs(phase3(k-1))) + phase_accu3(k-1);
   phase_accu3_0(k) = abs(abs(phase3_0(k)) - abs(phase3_0(k-1))) + phase_accu3_0(k-1);
end

%Calculate the effective index
n_eff1(k) = (c/L)*(phase_accu1(k)/(2*pi*fn(k)));
n_eff1_0(k) = (c/L)*(phase_accu1_0(k)/(2*pi*fn(k)));

n_eff2(k) = (c/L)*(phase_accu2(k)/(2*pi*fn(k)));
n_eff2_0(k) = (c/L)*(phase_accu2_0(k)/(2*pi*fn(k)));

n_eff3(k) = (c/L)*(phase_accu3(k)/(2*pi*fn(k)));
n_eff3_0(k) = (c/L)*(phase_accu3_0(k)/(2*pi*fn(k)));
end

%Calculate the gain profile and then nKK
for h = 1:length(frequency)
fn = frequency(h);
g_f(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - fB * 1e12)^2 + (FWHM_f/2)^2);
n_im(h) = (-1/2) * g_f(h) * c / (2 * pi * fn * 1e12);
delta_n_re(h) = - alpha_h * (fn - fB) * 1e12 / FWHM_f * n_im(h);
n_kk(h) = n_offset + delta_n_re(h);
end



%Plot figures
% figure (1);
% plot(frequency,n_avg*(n_eff1./n_eff1_0)); 
% xlim([2.5 3.3]);
% ylim([3.55 3.66]);
% xlabel('frequency THz');                                                   
% ylabel('Effective index');
% title('Periodic');
% 
% figure (2);
% plot(frequency,n_avg*(n_eff2./n_eff2_0)); 
% xlim([2.5 3.3]);
% ylim([3.55 3.66]);
% xlabel('frequency (THz)');                                                   
% ylabel('Effective index');
% title('Single-defect');
% 
% figure (3);
% plot(frequency,n_avg*(n_eff3./n_eff3_0)); 
% xlim([2.5 3.3]);
% ylim([3.55 3.66]);
% xlabel('frequency THz');                                                   
% ylabel('Effective index');
% title('Aperiodic');
% 
% figure (4);
% plot(frequency,g_f,'k');
% xlim([2.5 3.3]);
% ylim([0 1.1]);
% yticks(0:0.1:1.1);
% xlabel('frequency THz');
% ylabel('Gain (gLc)');
% title('Lorentzian gain profile');
% 
% figure (5);
% plot(frequency,n_kk,'r');
% xlim([2.5 3.3]);
% ylim([3.6025 3.607]);
% yticks(3.6025:0.0005:3.607);
% xlabel('frequency THz');
% ylabel('Refractive index (nkk)');
% title('Lorentzian gain profile');

figure (6);
yyaxis left; %Plot data on the left Y axis
plot(frequency,g_f*Lc,'k');
ylim([0 1.1]);
yticks(0:0.1:1.1);
ylabel('Gain (gLc)');

yyaxis right; %Plot data on the right Y axis
plot(frequency,n_kk,'r');
ylim([3.6025 3.607]);
yticks(3.6025:0.0005:3.607);
ylabel('Refractive index (nkk)');

%Set the x axis
xlim([2.5 3.3]);
xlabel('Frequency THz');
title('Lorentzian gain profile');
