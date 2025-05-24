clear all;

% PARAMETERS
n = [1, 1, 1, 1, 1.732, 3];     % n_air | n2 | n3 | ... | nN | n_substrate
d = [79.4, 79.4, 79.4, 79.4];             % thicknesses of layers (only for internal layers, no d for first and last)

lambda_range = 30:1:1400;       % Wavelength range in nm
num_layers = length(d);

% Phase shifts δ_j = 2π * n_j * d_j / λ for each layer and wavelength
delta = zeros(num_layers, length(lambda_range));
for j = 1:num_layers
    delta(j,:) = 2 * pi * n(j+1) * d(j) ./ lambda_range;  % n(j+1) is the j-th physical layer
end

r_total = zeros(size(lambda_range));
t_total = zeros(size(lambda_range));

% Loop over wavelengths
for idx = 1:length(lambda_range)
    lambda = lambda_range(idx);

    % Start from identity matrix
    T = eye(2);

    % Loop over each layer and interface
    for j = 1:num_layers
        % Fresnel interface j → j+1
        nj = n(j);
        nj1 = n(j+1);
        r = (nj - nj1) / (nj + nj1);
        t = 2 * nj / (nj + nj1);
        M_interface = (1 / t) * [1, r; r, 1];

        % Phase matrix for layer j (after interface j)
        phi = delta(j, idx);
        P = [exp(-1i * phi), 0;
             0, exp(1i * phi)];

        % Multiply T matrix: T = T * M_interface * P
        T = T * M_interface * P;
    end

    % Last interface (layer N to substrate)
    n_last = n(end-1);
    n_substrate = n(end);
    r_last = (n_last - n_substrate) / (n_last + n_substrate);
    t_last = 2 * n_last / (n_last + n_substrate);
    M_last = (1 / t_last) * [1, r_last; r_last, 1];
    T = T * M_last;

    % Calculate total r, t
    r_total(idx) = T(2,1) / T(1,1);
    t_total(idx) = 1 / T(1,1);
end

% INTENSITIES
I_r = abs(r_total).^2;
I_t = abs(t_total).^2 .* (n(end)/n(1));  % Flux correction for final substrate and initial medium

% PLOTTING
figure('Color', 'white');
plot(lambda_range, I_t, 'Color', [0.9 0.3 0.6], 'LineWidth', 2, 'DisplayName', 'Transmission'); % blue
hold on;
plot(lambda_range, I_r, 'Color', [0.5 0.2 0.7], 'LineWidth', 2, 'DisplayName', 'Reflection');  % red

xlabel('Wavelength (nm)', 'FontSize', 12);
ylabel('Normalized Intensity', 'FontSize', 12);
title('Transmission and Reflection vs. Wavelength', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 12);

set(gca, 'FontSize', 14, 'Box', 'on');
grid on;
