% Load the data from the file
data = load('output_data.txt');

% Extract the value 
sine = data(:, 1);  % Extract all rows in the first columns
cosine = data(:, 2);    %  Extract all rows in the second columns

sine = sine - mean(sine);
cosine = cosine - mean(cosine);

% Combine them into complex valued samples - sin to the real part, cosine
% to the imaginary part
complex_value = sine + 1j * cosine;

N = length(complex_value);


% Perform FFT
fft_output = fft(complex_value);



Fs = 100e6;
frequency_axis = Fs * (0:N-1) / N / 1e6;

figure;
plot(frequency_axis(1:N/2), abs(fft_output(1:N/2)), 'k', 'LineWidth', 1);
xlim([0, 5]);  % Focus on the range 0 to 5 MHz
xlabel('Frequency(MHz)');
ylabel('Absolute value');
title('FFT output');


