clear all;
fs = 0; fe = 2; 
Np = 1e4;
fc = 2.9;
frequency = fs*fc:(fe-fs)*fc/Np:fe*fc;

%The length of the whole structure
N = 203;
lamda = 14.3e-6;
Lg = N * lamda;
c = 3e8;

%Sets the refractive indices
kLg = 8;
n_eff = 3.605;     
delta_n = kLg / N * n_eff / 2;
n1 = n_eff + delta_n;
n2 = n_eff - delta_n;

%Calculate the interface constants S
S12=(1/(2*n1))*[n1+n2 n1-n2;n1-n2 n1+n2];
S21=(1/(2*n2))*[n2+n1 n2-n1;n2-n1 n2+n1];

% Set the gain coefficient g 
g = 2 / Lg; % gLg = 2

for i = 1:length(frequency)   

f = 1e12 * frequency(i);    

% Calculate the optical thickness
beta_AL_1 = 2 * pi * (f / c) * n_eff - 1i * g / 2;
beta_AL_2 = 2 * pi * (f / c) * n_eff - 1i * g / 2;

theta1 = beta_AL_1 * (lamda / 2);     %For slot d = lamda / 2  
theta2 = beta_AL_2 * (lamda / 2);     %For slot d = lamda / 2   
theta2_2 = beta_AL_2 * (lamda);     %For slot d = 1.5 * lamda   
theta2_3 = beta_AL_2 * (2 * lamda);     %For slot d = 2.5 * lamda   

%Calculate the medium constant P
Pn1 = [exp(-1i*theta1) 0;0 exp(1i*theta1)];
Pn2 = [exp(-1i*theta2) 0;0 exp(1i*theta2)];
Pn2_2 = [exp(-1i*theta2_2) 0;0 exp(1i*theta2_2)];
Pn2_3 = [exp(-1i*theta2_3) 0;0 exp(1i*theta2_3)];

%Construct the M matrix of the whole strcture
M = (Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^6 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^11 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^8 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^12*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^7*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_3*S21)*Pn1;

%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);

%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance(i)=Rpow;
end

% Set the gain coefficient g 
g = 0; % gLg = 0

for i = 1:length(frequency)   

f = 1e12 * frequency(i);    

% Calculate the optical thickness
beta_AL_1 = 2 * pi * (f / c) * n_eff - 1i * g / 2;
beta_AL_2 = 2 * pi * (f / c) * n_eff - 1i * g / 2;

theta1 = beta_AL_1 * (lamda / 2);     %For slot d = lamda / 2  
theta2 = beta_AL_2 * (lamda / 2);     %For slot d = lamda / 2   
theta2_2 = beta_AL_2 * (lamda);     %For slot d = 1.5 * lamda   
theta2_3 = beta_AL_2 * (2 * lamda);     %For slot d = 2.5 * lamda   

%Calculate the medium constant P
Pn1 = [exp(-1i*theta1) 0;0 exp(1i*theta1)];
Pn2 = [exp(-1i*theta2) 0;0 exp(1i*theta2)];
Pn2_2 = [exp(-1i*theta2_2) 0;0 exp(1i*theta2_2)];
Pn2_3 = [exp(-1i*theta2_3) 0;0 exp(1i*theta2_3)];

%Construct the M matrix of the whole strcture
M = (Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^6 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^11 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^2*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^4 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^2 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^8 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^12*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)*(Pn1*S12*Pn2_2*S21)^2*(Pn1*S12*Pn2*S21)^5 ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^7*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21) ...
    *(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^5*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2*S21)^3 ...
    *(Pn1*S12*Pn2_3*S21)*Pn1;

%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);

%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance_0(i)=Rpow;
end


figure;
plot(frequency, Reflectance, 'r', 'LineWidth', 1);
hold on;
plot(frequency, Reflectance_0, 'b', 'LineWidth', 1);
xlim([2.76 3.06]);
set(gca, 'YScale', 'log');
xlabel('frequency THz');                                                   
ylabel('Power reflectivity');






