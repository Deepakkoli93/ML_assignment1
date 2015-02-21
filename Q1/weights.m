function w = weights(datapoints, querypoint, t)

m = length(datapoints);

w = zeros(m, m);

for i=1:m
 weight = exp(-(((querypoint - datapoints(i))^2)/(2*t*t)));
 w(i,i) = weight;
 end;

