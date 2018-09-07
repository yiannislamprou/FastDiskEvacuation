% Fast Disk Evacuation: New Both-Explore Lower Bound Calculations
% Ioannis Lamprou
% March 2015

 % flag = 0;

for s = 1.0:0.1:7.2,
	%fprintf('s = %f\n', s);

	y_min = max((pi-s+1)/(s+1), 0); % y lower bound 
	y_max = (2*pi-s+1)/(s+1); % y upper bound - whole boundary can be explored in time 1 + y_max

	max_y = -1;
	argmax_y = [-1, -1];

	for y = y_min:0.001:y_max,
		min_k = 1000.0;
		argmin_k = -1;
		for k = 0:0.001:y,
			y_min_k = max(0, (pi-s+1+k)/(s+1));
			if y >= y_min_k,
				fast = 2*sin((s-1+(s+1)*y - k)/2)/s;
				slow = sqrt(sin((s-1+(s+1)*y - k)/2)^2 + max(1 - abs(cos((s-1+(s+1)*y - k)/2)) - k, 0)^2);
				% explore = (2*pi - s + 1 - (s+1)*y + k)/(s+1);
				%d = max(max(fast, slow), explore);
				d = max(fast, slow);
				if 1 + y + d < min_k,
					min_k = 1 + y + d;
					argmin_k = k;
				end
			end	
		end
		if min_k > max_y && min_k < 1000.0,
			max_y = min_k;
			argmax_y = [y, argmin_k];
		end
		explore = (2*pi - s + 1 - (s+1)*y)/(s+1);
		if 1 + y + explore > max_y, 
			max_y = 1+ y + explore;
			argmax_y = [y, -1];
		end	
	end
	
	fprintf('(%1.2f, %1.3f)', s, max_y);
	%fprintf('s = %1.2f, max = %1.4f, y = %1.4f, k = %1.4f\n', s, max_y, argmax_y(1), argmax_y(2));

	% if ~flag && s >= 2 && max_y > sqrt(1-4/s^2) + (1 + 2*acos(-2/s))/s, % gives 2.75
		% fprintf('First time OVER is %f\n', s);
		% flag = 1;
	% end
	
	% if s < 2 && max_y - 2*sqrt(1-s^2/(s+1)^2)/s - (-s + 2*acos(-s/(s+1))+1)/(s+1) - 1 > 10^(-10),
		% fprintf('Greater than old LB: s = %f, old = %f, new = %f\n', s, 2*sqrt(1-s^2/(s+1)^2)/s + (-s + 2*acos(-s/(s+1))+1)/(s+1) + 1, max_y);
	% end	

end