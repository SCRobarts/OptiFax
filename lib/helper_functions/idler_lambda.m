function l_i = idler_lambda(l_p,l_s)

c = 299792458;				% Speed of light (m/s)

f_p = c ./ l_p;
f_s = c ./ l_s;
f_i = f_p - f_s;

l_i = c ./ f_i;

end