function I = normDB(I,Imax,DBfloor)
	arguments
		I
		Imax = I;
		DBfloor = -60;
	end
		
	I	= I./ max(Imax(1,:));
		
	I	= pow2db(abs(I));
		
	I(I < DBfloor) = DBfloor;

	end