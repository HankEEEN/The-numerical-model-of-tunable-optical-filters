clear all; 

delta_f0 = (2.85 - 2.8) / 14;
[Reflectance1] = ReCal(2.85);
[Reflectance2] = ReCal(2.85 - 1 * delta_f0);
[Reflectance3] = ReCal(2.85 - 2 * delta_f0);
[Reflectance4] = ReCal(2.85 - 3 * delta_f0);
[Reflectance5] = ReCal(2.85 - 4 * delta_f0);
[Reflectance6] = ReCal(2.85 - 5 * delta_f0);
[Reflectance7] = ReCal(2.85 - 6 * delta_f0);
[Reflectance8] = ReCal(2.85 - 7 * delta_f0);
[Reflectance9] = ReCal(2.85 - 8 * delta_f0);
[Reflectance10] = ReCal(2.85 - 9 * delta_f0);
[Reflectance11] = ReCal(2.85 - 10 * delta_f0);
[Reflectance12] = ReCal(2.85 - 11 * delta_f0);
[Reflectance13] = ReCal(2.85 - 12 * delta_f0);
[Reflectance14] = ReCal(2.85 - 13 * delta_f0);
[Reflectance15] = ReCal(2.85 - 14 * delta_f0);




figure;
plot_figure_red(Reflectance1);
plot_figure_black(Reflectance2);
plot_figure_black(Reflectance3);
plot_figure_black(Reflectance4);
plot_figure_black(Reflectance5);
plot_figure_black(Reflectance6);
plot_figure_black(Reflectance7);
plot_figure_black(Reflectance8);
plot_figure_black(Reflectance9);
plot_figure_black(Reflectance10);
plot_figure_black(Reflectance11);
plot_figure_black(Reflectance12);
plot_figure_black(Reflectance13);
plot_figure_black(Reflectance14);
plot_figure_black(Reflectance15);


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
% Centre frequency of the grating filter in THz
fc = 2.9097;
% Frequency vector
frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;

% Refractive indices initialisation
% Set the kappa
kLg = 8;
% Set the effective index
n_eff = 3.605;
% Set the refractive indices according to kLg - Carroll delta_n definition
% P41
delta_n = kLg / N * n_eff / 2;
n1 = n_eff + delta_n;
n2 = n_eff - delta_n;


% Gain profile and complex refractive index nkk calculation
% Initialise all the parameters
g0 = 4 / Lc;
alpha_h = 4;

% With the gain dispersion - the modal gain follows the Lorentzian
% function, and the effective refractive index varies

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
beta_1 = 2 * pi * (fn / c) * n_kk(h) - 1i * g_f(h) / 2;
beta_2 = 2 * pi * (fn / c) * n_kk(h) - 1i * g_f(h) / 2;

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
xlim([2.838 2.846]);
set(gca, 'YScale', 'log');
% ylim([10^-0.5 10^5])
xlabel('Frequency THz');                                                   
ylabel('Intensity');

end

function [] = plot_figure_red (Reflectance)

Np = 1e5; %Sampling points

% Frequency Initialisation
fs = 0; fe = 2; %Set the range of frequency
fc = 2.9097; % Centre frequency in THz

frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;
plot(frequency, Reflectance, 'r', 'LineWidth', 1);
hold on;
xlim([2.838 2.846]);
set(gca, 'YScale', 'log');
% ylim([10^-0.5 10^5])
xlabel('Frequency THz');                                                   
ylabel('Intensity');

end





