load hw1.mat;

% % Prompt the user to input values for X, Y, and Z
% X = input('Enter the value of X: ');
% Y = input('Enter the value of Y: ');
% Z = input('Enter the value of Z: ');

% Value Initialization
best_X = 0;
best_Y = 0;
best_Z = 0;
best_snr = -Inf;

% Range for X, Y, Z
X_range = 1:10;
Y_range = 1:6;
Z_range = 1:9;

for X = X_range
    for Y = Y_range
        for Z = Z_range

            % Convert input signal and FIR taps into fixed-point
            fixed_signal = float_to_fixed(my_signal, 11, X);
            fixed_fir = float_to_fixed(my_fir_h, 7, Y);
            
            % floating point convolution
            result_float = conv(my_signal, my_fir_h);
            
            % fixed point convolution
            result_fixed = zeros(1, length(fixed_signal));
            
            for n = 1:length(fixed_signal)
                for k = 1:length(fixed_fir)
                    if n-k+1 > 0
                        mult_fixed = fixed_mul(fixed_signal(n-k+1), 11, X, fixed_fir(k), 7, Y, 10, Z);
                        result_fixed(n) = fixed_add(result_fixed(n), 10, Z, mult_fixed, 10, Z, 10, Z);
                    end
                end
            end
            
            % convert fixed-point output to the floating point output
            float_from_fixed = result_fixed * 2 ^ (-Z);
            
            % Plot the histogram of the error
            error = result_float(1:length(float_from_fixed)) - float_from_fixed;

            % Report the output SNR - comparing the same portion of the floating-point
            % and fixed-point signals, thus no need to normalize using mean. Sum will
            % be used here instead
            power_signal = sum(result_float .^ 2);
            power_error = sum(error .^ 2);
            snr = 10 * log10(power_signal / power_error);
%             fprintf('Output SNR: %.4f dB\n', snr);
            
            % Pick the best SNR
            if snr > best_snr
                best_X = X;
                best_Y = Y;
                best_Z = Z;
                best_snr = snr;
            end


        end
    end
end


fprintf('Output best_X: %d \n', best_X);
fprintf('Output best_Y: %d \n', best_Y);
fprintf('Output best_Z: %d \n', best_Z);
fprintf('Output SNR: %.4f dB\n', best_snr);


% Once we get the optimal values of X, Y, Z
X = best_X;
Y = best_Y;
Z = best_Z;

            % Convert input signal and FIR taps into fixed-point
            fixed_signal = float_to_fixed(my_signal, 11, X);
            fixed_fir = float_to_fixed(my_fir_h, 7, Y);
            
            % floating point convolution
            result_float = conv(my_signal, my_fir_h);
            
            % fixed point convolution
            result_fixed = zeros(1, length(fixed_signal));
            
            for n = 1:length(fixed_signal)
                for k = 1:length(fixed_fir)
                    if n-k+1 > 0
                        mult_fixed = fixed_mul(fixed_signal(n-k+1), 11, X, fixed_fir(k), 7, Y, 10, Z);
                        result_fixed(n) = fixed_add(result_fixed(n), 10, Z, mult_fixed, 10, Z, 10, Z);
                    end
                end
            end
            
            % convert fixed-point output to the floating point output
            float_from_fixed = result_fixed * 2 ^ (-Z);
            
            % Plot the histogram of the error
            error = result_float(1:length(float_from_fixed)) - float_from_fixed;

            % Report the output SNR - comparing the same portion of the floating-point
            % and fixed-point signals, thus no need to normalize using mean. Sum will
            % be used here instead
            power_signal = sum(result_float .^ 2);
            power_error = sum(error .^ 2);
            snr = 10 * log10(power_signal / power_error);
            fprintf('Output SNR: %.4f dB\n', snr);


% Plot that compares fixed point and floating point FIR results
figure(1);
plot(1:length(float_from_fixed), float_from_fixed, 'k', 'LineWidth', 1);
hold on;
plot(1:length(result_float), result_float, 'r', 'LineWidth', 1);
hold off;

legend('float from fixed', 'float');
xlabel('Sample');
ylabel('Value');
title('Comparison');


figure(2);
histogram(error);
xlabel('Error');
ylabel('Frequency (The number of samples)');
title('Difference between the floating point output and fixed point output');


