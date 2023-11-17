function t_field = k2t(field)

N = length(field);
% t_field = N*(fftshift(ifft(field)));
% t_field = (fftshift(ifft(field)));
t_field = N*(ifft(ifftshift(field)));

end