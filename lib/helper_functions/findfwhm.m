function y=findfwhm(x,q)
%finds fwhm of a well behaved function q on an abscissa of x
q=q/max(q);
p=find(q>0.5);
p0=min(p)-1;
p1=min(p);
p2=max(p);
p3=max(p)+1;

q0=q(p0);
q1=q(p1);
q2=q(p2);
q3=q(p3);

s01=(q1-q0)/(p1-p0);
s23=(q3-q2)/(p3-p2);

p01=p0+(0.5-q0)/s01;
p23=p2+(0.5-q2)/s23;

w=interp1(1:length(x),x,[p01 p23]);
if length(w) < 2
	y = range(x);
else
	y=abs(w(2)-w(1));
end
end
