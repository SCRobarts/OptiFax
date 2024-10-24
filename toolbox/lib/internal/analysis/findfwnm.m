function fwnm=findfwnm(x,y,n)
sz = size(y);
if sz(1) > sz(2)
	y = y.';
end
%finds full width nth maximum of a well behaved function q on an abscissa of x
%i.e. width over which q > max(q)*n
fwnm = zeros(min(sz),1);
for ii = 1:min(sz)
	q = y(ii,:);
	q=q./max(q);
	p=find(q>n);
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
	
	% p01=p0+(0.5-q0)/s01;
	% p23=p2+(0.5-q2)/s23;

	p01=p0+(n-q0)/s01;
	p23=p2+(n-q2)/s23;
	
	w=interp1(1:length(x),x,[p01 p23]);
	if length(w) < 2
		fwnm(ii) = range(x);
	else
		fwnm(ii)=abs(w(2)-w(1));
	end
end
end
