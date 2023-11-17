function k_field = t2k(field)

N = length(field);
% k_field = fft(ifftshift(field))./(N/1);
% k_field = fft(ifftshift(field));
k_field = fftshift(fft(field))./(N/1);

end