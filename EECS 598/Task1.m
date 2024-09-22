% Load the data
load hw1.mat;

% Determine the max and min values of the signal vector
max_signal = max(my_signal);
min_signal = min(my_signal);

% Calculate the uniform quantization step
delta = (max_signal - min_signal) / 32;

% Signal quantization
quantized_signal = min_signal + round((my_signal - min_signal) / delta) * delta;

% SQNR calculation
quantized_error = my_signal - quantized_signal;

% Since we are dealing with discrete signals, unlike the integration process in the
% definition of SQNR, a summation method will take the place of continuous
% integration
signal_p = mean(my_signal.^2);
noise_p = mean(quantized_error.^2);
sqnr = 10 * log10(signal_p / noise_p);

% Plot the original and quantized signal
figure(1);
plot(my_signal, 'k', 'LineWidth', 1);
hold on;
plot(quantized_signal, 'r', 'LineWidth', 1);
hold off;

% legend, label and title
legend('Original Signal', 'Quantized Signal');
xlabel('Sample');
ylabel('Value');
title('Comparison');

% Plot the PDF of quantization error
figure(2);
histogram(quantized_error, 'Normalization', 'pdf');
xlabel('Quantization Error');
ylabel('Probability Density');
title('PDF of Quantization Error');

% Report the corresponding 32 quantization levels and the best SQNR 
quantized_levels = min_signal:delta:max_signal; 
disp(quantized_levels);
disp(['Best SQNR achieved: ', num2str(sqnr), ' dB']);
