function [clamped_vals] = my_val_clamp(vals, smoothness)
    z = ((vals./ smoothness) + 1) ./ 2;
    clamped_vals = 1 * (z >= 1) + 0 * (z <= 0) + ...
                    z .* ((z < 1) & (z > 0));
