% Creating a function that converts a floating point number to a
% fixed-point quantized value

function fixed_mantissa = float_to_fixed(in, N, R)

% scale up using a factor
    scaleup = in * 2 ^ R;

% round to the nearest integer
    fixed_mantissa = round(scaleup);
    
% Check for appropriate saturation
    maximum = 2 ^ (N - 1) - 1;
    minimum = -2 ^ (N - 1);
    fixed_mantissa = max(minimum, min(maximum, fixed_mantissa));

end












