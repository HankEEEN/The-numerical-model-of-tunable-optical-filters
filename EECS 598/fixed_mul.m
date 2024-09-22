% Creating a function for fixed-point arithmetic operations -
% multiplication
function fixed_multiplication = fixed_mul(mantissa1, N1, R1, mantissa2, N2, R2, N3, R3)
    match_result = mantissa1 * mantissa2;
    result = match_result * 2 ^ (R1 + R2 - R3);

% Truncate the result 
    truncate_result = round(result);

% Saturation check
    maximum_value = 2 ^ (N3 - 1) - 1;
    minimum_value = -2 ^ (N3 - 1);

% Output the result
    fixed_multiplication = max(minimum_value, min(maximum_value, truncate_result));

end 
