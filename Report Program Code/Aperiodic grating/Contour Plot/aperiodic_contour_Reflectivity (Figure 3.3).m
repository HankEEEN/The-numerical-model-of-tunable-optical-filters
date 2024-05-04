clear all;

fs = 0; fe = 2; 
gs = 0; ge = 6;
Np = 3000;
Normalized_frequency = fs:(fe-fs)/Np:fe;
gL = gs:(ge-gs)/Np:ge;

% Speed of light
c = 3e8;
%Sets the refractive indices
n_eff = 3.644;
n1 = 3.69866; n2 = 3.58934; %kL = 1.5, neff = 3.644

% Set the central frequency of the grating structure fB
fB = 4.325; % THz

% Structure length
lamda = 9.51e-6;
L = 50 * lamda;

%Calculate the interface constants S
S12= (1 / (2 * n2)) * [n2 + n1 n2 - n1 ; n2 - n1 n2 + n1];
S21= (1 / (2 * n1)) * [n1 + n2 n1 - n2 ; n1 - n2 n1 + n2];

for j = 1:length(gL)
%Set the gain coefficient g
g = gL(j) / L;
for i = 1:length(Normalized_frequency)                                              
fn = fB * Normalized_frequency(i) * 1e12; 
% Calculate the propagation constant
beta = 2 * pi * (fn / c) * n_eff - 1i * g / 2;
%Calculate the optical thickness
theta1 = beta * (lamda / 2);    %For slot  lamda/2
theta2_1 = beta * (lamda / 2);    %For lamda lamda/2
theta2_2 = beta * (3 * lamda / 2);    %For 2*lamda 3lamda/2
theta2_3 = beta * (2 * lamda);    %For 2.5*lamda 2lamda
theta2_4 = beta * (3 * lamda);    %For 3.5*lamda 3lamda

%Calculate the medium constant P
Pn1 = [exp(-1i*theta1) 0;0 exp(1i*theta1)];
Pn2_1 = [exp(-1i*theta2_1) 0;0 exp(1i*theta2_1)];
Pn2_2 = [exp(-1i*theta2_2) 0;0 exp(1i*theta2_2)];
Pn2_3 = [exp(-1i*theta2_3) 0;0 exp(1i*theta2_3)];
Pn2_4 = [exp(-1i*theta2_4) 0;0 exp(1i*theta2_4)];

%Construct the M matrix of the whole strcture
M_left = ((Pn1*S12*Pn2_1*S21)^13);
M_mid = (Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2_1*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2_4*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2_1*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2_1*S21)*(Pn1*S12*Pn2_1*S21)*(Pn1*S12*Pn2_1*S21)*(Pn1*S12*Pn2_3*S21)*(Pn1*S12*Pn2_2*S21)*(Pn1*S12*Pn2_2*S21);
M_right = ((Pn1*S12*Pn2_1*S21)^14);
M = M_left * M_mid * M_right;

%Calculate the reflection coefficient
Ramp = M(2,1) / M(1,1);
%Calculate the power reflection coefficient
Rpow = (abs(Ramp)).^2;
Reflectance(i)=Rpow;
Z(j,i) = Reflectance(i);
end
end



%Create a contour plot
figure;
contour(Normalized_frequency, gL, Z, [0.01, 0.1, 1, 10, 100, 1000, 10000], 'k', 'LineWidth', 1, 'ShowText', 'on');


%Set the limits of the axes
xlim([0.92 1.08]);
ylim([0 6]);

%Lable and title the plot
xlabel('Normalized frequency');
ylabel('Gain (gL)');
title('Contour Plot of Reflectivity');
