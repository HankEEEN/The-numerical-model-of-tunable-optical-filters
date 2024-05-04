clear all; 

[Reflectance1] = ReCal(2.89);
[Reflectance2] = ReCal(2.803);

figure;
plot_figure_red(Reflectance1);
plot_figure_black(Reflectance2);


% Function for calculating the reflectance
function [Reflectance] = ReCal (f0)

% Sampling points
Np = 1e5;

% Speed of light
c = 3e8;

% Structure parameters - AL
N = 203;
lamda = 14.3e-6;
Lg = N * lamda;

% Structure parameters - AL + FP facets
Lc = 5e-3;

% Frequency tuning initialisation
% Set the range of frequency
fs = 0; fe = 2; 
% Centre frequency in THz
fc = 2.9097;
% Frequency vector
frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;

% Refractive indices initialisation
% Set the kappa
kLg = 8;
% Set the effective index
n_eff = 3.605;
% Set the refractive indices according to kLg
delta_n = kLg / N * n_eff / 2;
n1 = n_eff + delta_n;
n2 = n_eff - delta_n;


% Gain profile and complex refractive index nkk calculation
% Initialise all the parameters
g0 = 4 / Lc;
alpha_h = 4;

% Without the gain dispersion - the modal gain follows the Lorentzian
% function, but the effective refractive index 3.605 keeps constant

%Set the FWHM - full width at half maximum
FWHM_f = 300e9;

for h = 1:length(frequency)
fn = frequency(h);
g_f(h) = g0 * (FWHM_f/2)^2 / ((fn * 1e12 - f0 * 1e12)^2 + (FWHM_f/2)^2);
n_im(h) = (-1/2) * g_f(h) * c / (2 * pi * fn * 1e12);
delta_n_re(h) = - alpha_h * (fn - f0) * 1e12 / FWHM_f * n_im(h);
n_kk(h) = n_eff + delta_n_re(h);
end


% Calculate the interface constant S
S12 = (1 / (2 * n2)) * [n2 + n1 n2 - n1 ; n2 - n1 n2 + n1];
S21 = (1 / (2 * n1)) * [n1 + n2 n1 - n2 ; n1 - n2 n1 + n2];


% Calculate the propagation constant and optical thickness
for h = 1:length(frequency)                                              

fn = 1e12 * frequency(h) ;  
beta_1 = 2 * pi * (fn / c) * n_eff - 1i * g_f(h) / 2;
beta_2 = 2 * pi * (fn / c) * n_eff - 1i * g_f(h) / 2;

theta_1 = beta_1 * (lamda / 2) ;     %For slot d = lamda / 2  
theta_2 = beta_2 * (lamda / 2) ;     %For slot d = lamda / 2   
theta_2_2 = beta_2 * (lamda) ;     %For slot d = lamda   
theta_2_3 = beta_2 * (2 * lamda) ;     %For slot d = 2 * lamda  

% Calculate the medium constant P 
Pn1 = [exp(-1i*theta_1) 0;0 exp(1i*theta_1)];
Pn2 = [exp(-1i*theta_2) 0;0 exp(1i*theta_2)];
Pn3 = [exp(-1i*theta_2_2) 0;0 exp(1i*theta_2_2)];
Pn4 = [exp(-1i*theta_2_3) 0;0 exp(1i*theta_2_3)];

% Construct the M matrix of the whole strcture - B
M = (Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^6 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^11 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^4 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^8 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^12*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn3*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^7*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn3*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn4*S21)*Pn1;


% Calculate the reflection coefficient
Rre = M(2,1) / M(1,1);

% Calculate the power reflection coefficient
Rpow = (abs(Rre)).^2;
Reflectance(h) = Rpow;


end

end

% Function for figure plotting
function [] = plot_figure_black (Reflectance)

Np = 1e5; %Sampling points

% Frequency Initialisation
fs = 0; fe = 2; %Set the range of frequency
fc = 2.9097; % Centre frequency in THz

frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;
plot(frequency, Reflectance, 'k', 'LineWidth', 1);
hold on;
xlim([2.77 2.93]);
set(gca, 'YScale', 'log');
ylim([10^-0.5 10^5])
xlabel('Frequency THz');                                                   
ylabel('Reflectivity');

end

function [] = plot_figure_red (Reflectance)

Np = 1e5; %Sampling points

% Frequency Initialisation
fs = 0; fe = 2; %Set the range of frequency
fc = 2.9097; % Centre frequency in THz

frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;
plot(frequency, Reflectance, 'r', 'LineWidth', 1);
hold on;
xlim([2.77 2.93]);
set(gca, 'YScale', 'log');
ylim([10^-0.5 10^5])
xlabel('Frequency THz');                                                   
ylabel('Reflectivity');

end






