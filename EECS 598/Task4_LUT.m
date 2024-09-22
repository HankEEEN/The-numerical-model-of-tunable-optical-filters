% Define the number of entries
N = 256;

% Maximum positive value for saturation test
maxi = 127 / 2 ^ 7;

% LUT array
LUT = zeros(1, N);

for k = 0 : N-1
    theta = 2 * pi * k / N;
    cosine = cos(theta);

    % Saturation test
    if cosine > maxi
        cosine = maxi;
    end

    % Quantization to (8,7) fixed-point number
    if cosine >= 0 
        fixed_cos = round(cosine * 2 ^ 7);
    else
        % Obtain the 2's complement for negative cos values utilizing mod
        fixed_cos = mod(round(cosine * 2 ^ 7), 256);
    end

    % MATLAB characteristic - Index starts from 1
    LUT(k+1) = fixed_cos;
end


% Generate the output in preparation for the case statement in verilog
% program

fprintf('always @(*) begin \n');
fprintf('   case (entry)\n');
for k = 0 : N - 1
    fprintf('       8''d%d: data = 8''b%s;\n', k, dec2bin(LUT(k+1), 8));
end
fprintf('       default: data = 8''d0;\n');
fprintf('   endcase\n');
fprintf('end\n');

