% Fast Disk Evacuation: Fast-Chord Upper Bound
% Ioannis Lamprou
% March 2015

for s = 1.0:0.1:7.2,

	% Try different values based on x3
	min3 = 100;
	argmin3 = [-1, -1, -1, -1, -1, -1];
	
	for x3 = 0:0.01:(2*pi-s+1),%

		syms Sx1 Sx2 Sy;
		[v1, v2, v3] = vpasolve([Sx1 + Sy + x3 + s-1 == 2*pi, Sx2 == 2*sin((x3 + Sy)/2), Sx1 + Sx2 == s*Sy], [Sx1, Sx2, Sy]);
		
		x1 = double(v1);
		x2 = double(v2);
		y = double(v3);
		
		if x1 < 0 || x2 < 0 || y < 0,
			continue;
		end	
		
		% Compute upper bounds
	
		% Phase I: Before Slow explores
		maxI = -1;
		argmaxI = -1;
		for t = 1/s:0.01:1,
			if t + sqrt(1+t^2-2*t*cos(s*(1-t)+x1)) > maxI,
				maxI = t + sqrt(1+t^2-2*t*cos(s*(1-t)+x1));
				argmaxI = t;
			end	
		end 
		
		% Phase II: While Slow explores
		maxII = -1;
		argmaxII = -1;
		subcaseII = -1;
		
		% Phase IIa: While Fast in x1
		for t = 1:0.01:1 + x1/s,	
			dist = sqrt((cos(s*(t - 1/s)) - cos(s-1 + x1 + t - 1))^2 + (sin(s*(t - 1/s)) - sin(s-1 + x1 + t - 1))^2);
			if t + dist > maxII,
				maxII = t + dist;
				argmaxII = t;
				subcaseII = 1;
			end	
		end
		
		% Phase IIb: While Fast in x2
		for t = 1 + x1/s:0.01:1 + x1/s + x2/s,	
			fast_x = cos(s-1 + x1) + s*(t - 1 - x1/s)*(1 - cos(s-1 + x1))/x2;
			fast_y = sin(s-1 + x1) + s*(t - 1 - x1/s)*(-sin(s-1 + x1))/x2;
			dist = sqrt((fast_x - cos(s-1 + x1 + t - 1))^2 + (fast_y - sin(s-1 + x1 + t - 1))^2);
			if t + dist > maxII,
				maxII = t + dist;
				argmaxII = t;
				subcaseII = 2;
			end	
		end
		
		% Phase IIc: While Fast in x3
		for t = 1 + x1/s + x2/s:0.01:1 + x1/s + x2/s + x3/(s+1),	
			dist = sqrt((cos(2*pi - s*(t - 1 - x1/s - x2/s)) - cos(s-1 + x1 + t - 1))^2 + (sin(2*pi - s*(t - 1 - x1/s - x2/s)) - sin(s-1 + x1 + t - 1))^2);
			if t + dist > maxII,
				maxII = t + dist;
				argmaxII = t;
				subcaseII = 3;
			end	
		end
		
		% Max out of I, II
		if maxI > maxII,
			maximum = maxI;
			argmax = argmaxI;
			subcaseII = -1;
		else
			maximum = maxII;
			argmax = argmaxII;
		end	
		
		% Min over x3
		if maximum < min3,
			min3 = maximum;
			argmin3 = [argmax, x1, x2, x3, y, subcaseII];
		end	
		
	end
	
	fprintf('(%1.2f, %1.3f)', s, min3);
	%fprintf('Worst-case for s = %1.2f is %1.3f for t = %1.3f, x1 = %1.3f, x2 = %1.3f, x3 = %1.3f, y = %1.3f and subcase %d\n', s, min3, argmin3(1), argmin3(2), argmin3(3), argmin3(4), argmin3(5), argmin3(6));
	
end