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


% --------------------------- a --------------------------------
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

end;
disp('final values of theta and cost are as follows...');
theta
cost


% --------------------------------------------------------------
disp(' part a is done press enter to continue');
pause;
hold off;
% --------------------------- b --------------------------------
disp('plotting graphs ...');
pause(1);
plot(X(:,2),y,'o');
xlabel('area of the houses');
ylabel('prices of the houses');
hold;

% extracting two points to plot the line
hypo_x = [min(X(:,2))-0.1 ,max(X(:,2))+0.1 ];
hypo_y = theta(1) + theta(2)*hypo_x;
plot(hypo_x,hypo_y);

% ---------------------------------------------------------------
%disp(' part b is done press enter to continue');
pause;
hold off;
% --------------------------- c --------------------------------
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

mesh(theta0_vals, theta1_vals, J_vals);
rotate3d;
hold;
theta = zeros(size(X,2),1);


num_iter = 10000;
alpha = 0.01;
epsilon = 0.01;
cost_old=0;
cost=10000;
iter=0;
%repeat until convergence

while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	scatter3(theta(1),theta(2),cost,'x');
	iter=iter+1;
	view(-37.5,30 - iter/2);
pause (0.2);

end;


% --------------------------------------------------------------
disp(' part c is done press enter to continue');
pause;
hold off;
% --------------------------- d --------------------------------
%{
theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
%}
% --------------------------------------------------------------
%disp(' part d is done press enter to continue');
%pause;

hold off;
% --------------------------- d 0.1 --------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------

hold off;
% --------------------------- d 0.5 --------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------

hold off;
% --------------------------- d 0.9--------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------


hold off;
% --------------------------- d 1.3 --------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------


hold off;
% --------------------------- d 2.1 --------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------

hold off;
% --------------------------- d 2.5 --------------------------------
%{theta0_vals = linspace(-100, 100, 2000);
theta1_vals = linspace(-70, 70, 2000);
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
alpha = 2.1;
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
%}
% --------------------------------------------------------------











