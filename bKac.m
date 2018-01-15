function g = bKac(g,S,mask,partial)
    % block Kaczmarz method
    not_mask_index = find(mask==1);
    repeat = 100*ceil(1/partial);
    F_g = fftshift(fft2(g));
    theta_l = F_g./abs(F_g);
    total_delta = S.*mask.*theta_l-F_g;
    for i=1:repeat
    partial_index = randsample(not_mask_index,ceil(partial*numel(g)));
    partial_delta = zeros(size(g));
    partial_delta(partial_index) = total_delta(partial_index);
    delta_g = ifft2(ifftshift(partial_delta));
    g = g + 0.5*delta_g;
    end
    g = abs(g);
end