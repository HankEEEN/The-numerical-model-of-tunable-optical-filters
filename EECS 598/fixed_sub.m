% Creating a function for fixed-point arithmetic operations - subtraction
function fixed_subtraction = fixed_sub(mantissa1, N1, R1, mantissa2, N2, R2, N3, R3)
% Match the exponents of two operands with the output integer mantissa
    match_mantissa1 = mantissa1 * 2 ^ (R3 - R1);
    match_mantissa2 = mantissa2 * 2 ^ (R3 - R2);

    match_result = match_mantissa1 - match_mantissa2;

% Saturation check
    maximum_value = 2 ^ (N3 - 1) - 1;
    minimum_value = -2 ^ (N3 - 1);

% Output the result
    fixed_subtraction = max(minimum_value, min(maximum_value, match_result));


end
