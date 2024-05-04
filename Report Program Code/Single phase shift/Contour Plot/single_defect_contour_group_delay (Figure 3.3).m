clear all;

fs = 0; fe = 2; 
gs = 0; ge = 6;
Np = 3000;
Normalized_frequency = fs:(fe-fs)/Np:fe;
gL = gs:(ge-gs)/Np:ge;
delta_f = (fe-fs)/1e10;

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
theta2_2 = beta * (lamda);    %For 1.5*lamda lamda
%Calculate the medium constant P
Pn1 = [exp(-1i*theta1) 0;0 exp(1i*theta1)];
Pn2_1 = [exp(-1i*theta2_1) 0;0 exp(1i*theta2_1)];
Pn2_2 = [exp(-1i*theta2_2) 0;0 exp(1i*theta2_2)];
%M matrix of the uniform grating structure
M_left = ((Pn1*S12*Pn2_1*S21)^24);
M_mid = Pn1*S12*Pn2_2*S21;
M_right = ((Pn1*S12*Pn2_1*S21)^25);
M = M_left * M_mid * M_right;
%Calculate the transmission coefficient
Rtrans = 1 / M(1,1);

%Calculate the phase 
phase_1(j,i) = angle(Rtrans);

end
end

for j = 1:length(gL)
%Set the gain coefficient g
g = gL(j) / L;
for i = 1:length(Normalized_frequency)                                              
fn = fB * (Normalized_frequency(i) + delta_f) * 1e12; 
% Calculate the propagation constant
beta = 2 * pi * (fn / c) * n_eff - 1i * g / 2;
%Calculate the optical thickness
theta1 = beta * (lamda / 2);    %For slot  lamda/2
theta2_1 = beta * (lamda / 2);    %For lamda lamda/2
theta2_2 = beta * (lamda);    %For lamda lamda
%Calculate the medium constant P
Pn1 = [exp(-1i*theta1) 0;0 exp(1i*theta1)];
Pn2_1 = [exp(-1i*theta2_1) 0;0 exp(1i*theta2_1)];
Pn2_2 = [exp(-1i*theta2_2) 0;0 exp(1i*theta2_2)];
%M matrix of the uniform grating structure
M_left = ((Pn1*S12*Pn2_1*S21)^24);
M_mid = Pn1*S12*Pn2_2*S21;
M_right = ((Pn1*S12*Pn2_1*S21)^25);
M = M_left * M_mid * M_right;
%Calculate the transmission coefficient
Rtrans = 1 / M(1,1);

%Calculate the phase 
phase_2(j,i) = angle(Rtrans);

end
end

%Calculate the group delay using difference method
for r = 1:length(gL)
for k = 1:length(Normalized_frequency)
Z(r,k) = (phase_2(r,k) - phase_1(r,k)) / (2 * pi * (fB * delta_f * 1e12));
end
end



%Create a contour plot
figure;
contour(Normalized_frequency, gL, Z * 1e12, [1, 2.1, 4.6, 10, 21, 46, 100], 'k','ShowText','on');


%Set the limits of the axes
xlim([0.92 1.08]);
ylim([0 6]);

%Labels and title
xlabel('Normalized frequency');
ylabel('Gain (gL)');
title('Contour of Group Delay (ps)');




