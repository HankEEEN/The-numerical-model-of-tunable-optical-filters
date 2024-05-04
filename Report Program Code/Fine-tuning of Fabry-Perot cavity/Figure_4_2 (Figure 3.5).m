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
g0 = 1 / Lc;
alpha_h = 4;

%Frequency initialisation
%Set the range of frequency
fs = 0; fe = 2; 
%Set the FWHM - full width at half maximum
FWHM_f = 300e9;
%Centre frequency in THz
fB = 2.85;
fC = 2.95;

frequency_B = fs*fB:(fe-fs)*fB/Np:fe*fB;
frequency_C = fs*fC:(fe-fs)*fC/Np:fe*fC;

%Set the refractive indices
n_offset = 3.867;

%Calculate the gain profile and then nKK
for h = 1:length(frequency_B)
fn = frequency_B(h);
g_f_B(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - fB * 1e12)^2 + (FWHM_f/2)^2);
n_im_B(h) = (-1/2) * g_f_B(h) * c / (2 * pi * fn * 1e12);
delta_n_re_B(h) = - alpha_h * (fn - fB) * 1e12 / FWHM_f * n_im_B(h);
n_kk_B(h) = n_offset + delta_n_re_B(h);
end

for h = 1:length(frequency_C)
fn = frequency_C(h);
g_f_C(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - fC * 1e12)^2 + (FWHM_f/2)^2);
n_im_C(h) = (-1/2) * g_f_C(h) * c / (2 * pi * fn * 1e12);
delta_n_re_C(h) = - alpha_h * (fn - fC) * 1e12 / FWHM_f * n_im_C(h);
n_kk_C(h) = n_offset + delta_n_re_C(h);
end


figure (1);
plot(frequency_B, g_f_B * Lc, 'k', frequency_C, g_f_C * Lc, 'r');
title('Lorentzian gain profile');
xlim([2.75 3.05]);
xlabel('Frequency THz');
ylim([0 1.3]);
yticks(0.25:0.25:1.3);
ylabel('Gain (gLc)');


figure (2);
plot(frequency_B, n_kk_B, 'k', frequency_C, n_kk_C, 'r');
title('Calculated refractive indices');
xlim([2.75 3.05]);
xlabel('Frequency THz');
ylim([3.865 3.869]);
yticks(3.865:0.0005:3.869);
ylabel('Refractive index (nkk)');


