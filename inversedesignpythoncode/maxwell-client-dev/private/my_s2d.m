function [d] = s2d(s)
    if any(isinf(s))
        d = 1;
    else
        d = real(s(:));
    end
end


