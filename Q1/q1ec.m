X = load('q1x.dat');
y = load('q1y.dat');

m = length(y);

% normalizing
mu = mean(X(:,1));
sigma = std(X(:,1));
X = (X-mu)./sigma;

%adding the intercept term
X = [ones(m,1) X];

theta = zeros(size(X,2),1);

%----------------------------



% alpha = 0.1
disp('aplha = 0.1');
hold off;
theta0_vals = linspace(0, 12, 1000);
theta1_vals = linspace(0, 10, 1000);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];  
	  h = X*t;
e = h-y;
v_sqr = e.^2;
q=sum(v_sqr);
J=q/m;
J=J/2;  
	  J_vals(i,j) = J;
    end
end

contour(theta0_vals, theta1_vals, J_vals);

hold;
theta = zeros(size(X,2),1);


num_iter = 10000;
alpha = 0.1;
epsilon = 0.00001;
cost_old=0;
cost=10000;

%repeat until convergence

while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	scatter3(theta(1),theta(2),cost,'x');
pause (0.2);

end;

disp('press enter to continue');
pause;


% alpha = 0.5
disp('aplha = 0.5');

hold off;
theta0_vals = linspace(0, 12, 1000);
theta1_vals = linspace(0, 10, 1000);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];  
	  h = X*t;
e = h-y;
v_sqr = e.^2;
q=sum(v_sqr);
J=q/m;
J=J/2;  
	  J_vals(i,j) = J;
    end
end

contour(theta0_vals, theta1_vals, J_vals);

hold;
theta = zeros(size(X,2),1);


num_iter = 10000;
alpha = 0.5;
epsilon = 0.00001;
cost_old=0;
cost=10000;

%repeat until convergence

while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	scatter3(theta(1),theta(2),cost,'x');
pause (0.2);

end;


disp('press enter to continue');
pause;


% alpha = 0.9
disp('aplha = 0.9');

hold off;
theta0_vals = linspace(0, 12, 1000);
theta1_vals = linspace(0, 10, 1000);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];  
	  h = X*t;
e = h-y;
v_sqr = e.^2;
q=sum(v_sqr);
J=q/m;
J=J/2;  
	  J_vals(i,j) = J;
    end
end

contour(theta0_vals, theta1_vals, J_vals);

hold;
theta = zeros(size(X,2),1);


num_iter = 10000;
alpha = 0.9;
epsilon = 0.00001;
cost_old=0;
cost=10000;

%repeat until convergence

while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	scatter3(theta(1),theta(2),cost,'x');
pause (0.2);

end;


disp('press enter to continue');
pause;


% alpha = 1.3
disp('aplha = 1.3');

hold off;
theta0_vals = linspace(0, 12, 1000);
theta1_vals = linspace(0, 10, 1000);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
	  t = [theta0_vals(i); theta1_vals(j)];  
	  h = X*t;
e = h-y;
v_sqr = e.^2;
q=sum(v_sqr);
J=q/m;
J=J/2;  
	  J_vals(i,j) = J;
    end
end

contour(theta0_vals, theta1_vals, J_vals);

hold;
theta = zeros(size(X,2),1);


num_iter = 10000;
alpha = 1.3;
epsilon = 0.00001;
cost_old=0;
cost=10000;

%repeat until convergence

while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	scatter3(theta(1),theta(2),cost,'x');
pause (0.2);

end;
