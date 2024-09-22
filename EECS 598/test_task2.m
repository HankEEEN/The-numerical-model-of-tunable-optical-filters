% Convert floating point to fixed-point (example: 3.1416 in (8,4))
fixed_mantissa = float_to_fixed(32480, 8, 4);
disp(['Fixed-point mantissa: ', num2str(fixed_mantissa)]);

% Fixed-point addition example: (mantissa1: (8,4), mantissa2: (8,3)) => (8,4)
fixed_addition = fixed_add(50, 8, 4, 40, 8, 3, 8, 4);
disp(['Addition result: ', num2str(fixed_addition)]);

% Fixed-point subtraction example: (mantissa1: (8,4), mantissa2: (8,3)) => (8,4)
fixed_subtraction = fixed_sub(50, 8, 4, 40, 8, 3, 8, 4);
disp(['Subtraction result: ', num2str(fixed_subtraction)]);

% Fixed-point multiplication example: (mantissa1: (8,4), mantissa2: (8,3)) => (8,4)
fixed_multiplication = fixed_mul(50, 8, 4, 10, 8, 3, 8, 4);
disp(['Multiplication result: ', num2str(fixed_multiplication)]);