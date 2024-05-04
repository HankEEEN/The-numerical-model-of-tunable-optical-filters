clear all;

%Sampling points
Np = 1e6;

c = 3e8;

%Set the modal gain g0
Lc = 5e-3;
g0 = 1 / Lc;
alpha_h = 4;

%Frequency initialisation
%Set the range of frequency
fs = 0; fe = 2; 
%Set the FWHM - full width at half maximum
FWHM_f = 300e9;
%Centre frequency in THz
fc = 2.93;
f0_1 = 2.85;
f0_2 = 2.95;

frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;

%Set the refractive indices
n_eff = 3.867;
n1 = n_eff; n2 = n_eff; n0 = 1;


%Calculate the gain profile and then nKK
for h = 1:length(frequency)
fn = frequency(h);
g_f_1(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - f0_1 * 1e12)^2 + (FWHM_f/2)^2);
n_im_1(h) = (-1/2) * g_f_1(h) * c / (2 * pi * fn * 1e12);
delta_n_re_1(h) = - alpha_h * (fn - f0_1) * 1e12 / FWHM_f * n_im_1(h);
n_kk_1(h) = n_eff + delta_n_re_1(h);
end

for h = 1:length(frequency)
fn = frequency(h);
g_f_2(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - f0_2 * 1e12)^2 + (FWHM_f/2)^2);
n_im_2(h) = (-1/2) * g_f_2(h) * c / (2 * pi * fn * 1e12);
delta_n_re_2(h) = - alpha_h * (fn - f0_2) * 1e12 / FWHM_f * n_im_2(h);
n_kk_2(h) = n_eff + delta_n_re_2(h);
end




%Calculate the interface constants S
S01 = (1/(2*n0))*[n0+n1 n0-n1;n0-n1 n0+n1];
S10 = (1/(2*n1))*[n1+n0 n1-n0;n1-n0 n1+n0];
S02 = (1/(2*n0))*[n0+n2 n0-n2;n0-n2 n0+n2];
S20 = (1/(2*n2))*[n2+n0 n2-n0;n2-n0 n2+n0];

for h = 1:length(frequency)                                              

fn = 1e12 * frequency(h) ;  

%Calculate the optical thickness
beta_FP_1 = 2 * pi * (fn / c) * n_kk_1(h) - 1i * g_f_1(h) / 2;
theta_1 = beta_FP_1 * Lc  ;     

beta_FP_2 = 2 * pi * (fn / c) * n_kk_2(h) - 1i * g_f_2(h) / 2;
theta_2 = beta_FP_2 * Lc  ;         

%Calculate the medium constant P
Pn1_1 = [exp(-1i*theta_1) 0;0 exp(1i*theta_1)];

Pn1_2 = [exp(-1i*theta_2) 0;0 exp(1i*theta_2)];

%Construct the M matrix of the whole strcture
M_1 = S01 * Pn1_1 * S10;
M_2 = S02 * Pn1_2 * S20;

%Calculate the reflection coefficient
Ramp_1 = M_1(2,1) / M_1(1,1);
Ramp_2 = M_2(2,1) / M_2(1,1);
Rtrans_1 = 1 / M_1(1,1);
Rtrans_2 = 1 / M_2(1,1);

%Calculate the phase 
phase_1(h) = angle(Rtrans_1);
phase_2(h) = angle(Rtrans_2);

%Calculate the power reflection coefficient
Rpow_1 = (abs(Ramp_1)).^2;
Reflectance_1(h) = Rpow_1;

Rpow_2 = (abs(Ramp_2)).^2;
Reflectance_2(h) = Rpow_2;

end

for k = 2:length(frequency)
%When k=1, n_kk_1 and n_kk_2 are both equal to minus infinity, leading to
%phase(1) = NaN. Thus the case when k=1 should be ignored
if k==2
   phase_accu_1(k) = abs(phase_1(k));
   phase_accu_2(k) = abs(phase_2(k));
else
   phase_accu_1(k) = abs(abs(phase_1(k)) - abs(phase_1(k-1))) + phase_accu_1(k-1);
   phase_accu_2(k) = abs(abs(phase_2(k)) - abs(phase_2(k-1))) + phase_accu_2(k-1);
end

end


figure;
yyaxis left; %Plot data on the left Y axis
plot(frequency, Reflectance_1, 'k', 'LineWidth', 1);
hold on;
plot(frequency, Reflectance_2, 'r', 'LineWidth', 1);
hold off;
ylabel('Power reflectivity');

yyaxis right; %Plot data on the right Y axis
plot(frequency, phase_accu_1 / pi, 'k', 'LineWidth', 1);
hold on;
plot(frequency, phase_accu_2 / pi, 'r', 'LineWidth', 1); 
hold off;
ylabel('Transmission phase (pi)');


xlim([2.892 2.906]);
xlabel('Frequency THz');                                                   
title('Reflectance Spectrum of Fabry-Perot Cavity'); 



