X = load('q3x.dat');
y = load('q3y.dat');

m = length(y);

%adding the intercept term
X = [ones(m,1) X];


%--------------------- a -----------------------
theta = zeros(size(X,2),1);

% analytical solution of theta
theta = pinv(transpose(X)*X)*transpose(X)*y;
disp('theta is as follows');
theta
hold off;
plot(X(:,2),y,'.');
xlabel('x');
ylabel('y');
hold;

% extracting two points to define the line 
  plot_x = [min(X(:,2))-0.5,  max(X(:,2))+0.5];

% calculating decision boundary
  plot_y = theta(1) + theta(2)*plot_x;

plot(plot_x, plot_y);

%-----------------------------------------------
hold off;
disp(' press enter to continue...')
pause;

%--------------------- b -----------------------------------------------
disp('plotting graph');
pause(1);

datapoints = X(:,2);
querypoints = linspace(plot_x(1), plot_x(2), 100);

plot(X(:,2),y,'.');
xlabel('x');
ylabel('y');
hold;

t = 0.8;

for i=1:length(querypoints)
	weight = weights(datapoints,querypoints(i), t);
	theta = pinv(transpose(X) * weight * X) * transpose(X) * weight * y;
	prediction = transpose(theta) * [1;querypoints(i)];
	%pause (0.1);
	plot(querypoints(i),prediction,'x');
end

%-------------------------------------------------------------------------
hold off;
disp('press enter to continue to part c');
pause;
%--------------------- c -----------------------------------------------
taus = [0.1 0.3 2 10]; 
for z=1:4
	t = taus(z)
hold off;
figure();
disp('plotting graph');
pause(1);

datapoints = X(:,2);
querypoints = linspace(plot_x(1), plot_x(2), 100);

plot(X(:,2),y,'.');
xlabel('x');
ylabel('y');
title(t);
hold;

%t = 0.8;

for i=1:length(querypoints)
	weight = weights(datapoints,querypoints(i), t);
	theta = pinv(transpose(X) * weight * X) * transpose(X) * weight * y;
	prediction = transpose(theta) * [1;querypoints(i)];
	%pause (0.1);
	plot(querypoints(i),prediction,'x');
end
pause;
end;



