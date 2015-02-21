X = load('~/semester7/csl341/Assignment1/q1x.dat');
y = load('~/semester7/csl341/Assignment1/q1y.dat');

m = length(y);

% ------------------ 1 a --------------------

% normalizing
mu = mean(X(:,1));
sigma = std(X(:,1));
X = (X-mu)./sigma;


%adding the intercept term
X = [ones(m,1) X];

theta = zeros(size(X,2),1);

% X -> m * (n+1)
% y -> m * 1
% theta -> (n+1) * 1
hypothesis = X * theta;
error = sum(((y - hypothesis).^2)./2);
fname = '~/Desktop/sample_normalized.txt';
fid = fopen(fname,'w');

% perfect! giving optimum cost = 4.4770
num_iter = 10000;
alpha = 0.01;
epsilon = 0.00000001;
cost_old=100000;
cost=10000;
%repeat until convergence
%for i=1:num_iter
while (abs(cost_old-cost)>epsilon)
	cost_old=cost;
	hypothesis = X * theta;
	cost = (1/(2*m)) * sum(((hypothesis)-y).^2);
	theta = theta - (alpha * transpose(X) * (hypothesis - y))/m;
	fprintf(fid,'%s ',theta(1));
	fprintf(fid,'%s ',theta(2));
	fprintf(fid,'%s  \n',cost);
	scatter3(theta(1),theta(2),cost,'x');
pause (0.2);

end;
% ------------------ 1 a over --------------------

% ------------------ 1 b --------------------
%plotting graphs
plot(X(:,2),y,'o');
hold;
hypo_x = 4:2:24;
hypo_y = theta(1) + theta(2)*hypo_x;
plot(hypo_x,hypo_y);
% ---------------- 1 b over ------------------

% ------------------ 1 c --------------------
%plotting 3d mesh
theta0_vals = linspace(-3, 7, 2000);
theta1_vals = linspace(0, 1.2, 2000);
J_vals = zeros(length(theta0_vals), length(theta1_vals));

h = X*theta;
e = h-y;
v_sqr = e.^2;
q=sum(v_sqr);
J=q/m;
J=J/2;

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

min = 50;
for i = 1:2000
    for j = 1:2000

if(J_vals(i,j)<min)
min=J_vals(i,j);
end
end
end

%----------------------xxxxxxxxxxxxxx question 3 xxxxxxxx--------------------
X = load('q3x.dat');
y = load('q3y.dat');

m = length(y);

% ------------------ 3 a --------------------

%adding the intercept term
X = [ones(m,1) X];

theta = zeros(size(X,2),1);

theta = pinv(transpose(X)*X)*transpose(X)*y;

%plotting graphs
plot(X(:,2),y,'o');
hold;
hypo_x = -5:2:12;
hypo_y = theta(1) + theta(2)*hypo_x;
plot(hypo_x,hypo_y);
hold;
% ------------------ 3 a over --------------------

% ------------------ 3 b --------------------

datapoints = X(:,2);
querypoints = linspace(-5, 12, 100);

for i=1:length(querypoints)
weight = weights(datapoints,querypoints(i), 0.1);
theta = pinv(transpose(X) * weight * X) * transpose(X) * weight * y;
prediction = transpose(theta) * [1;querypoints(i)];
pause (0.1);
plot(querypoints(i),prediction,'x');
end

% ------------------ 3 b over --------------------

% ------------------- 2 a ----------------------

X = load('q2x.dat');
y = load('q2y.dat');

m = length(y);
%adding the intercept term
X = [ones(m,1) X];

theta = zeros(size(X,2),1);

n = length(theta);


oldtheta = 100.*ones(n,1);
epsilon = 0.01;
iter = 1;
while(norm(theta - oldtheta)>epsilon)
	iter = iter + 1
	oldtheta = theta;
	theta = theta - pinv(hessian(X,y,theta))*grad(X,y,theta);

   end; 
% ------------------- 2 a over ----------------------

% ------------------- 2 b ----------------------
pos = find(y==1); neg = find(y==0);

plot(X(pos,2),X(pos,3),'+');
hold
plot(X(neg,2),X(neg,3),'o');

% Only need 2 points to define a line, so choose two endpoints
    plot_x = [min(X(:,2))-0.1,  max(X(:,2))+0.1];

% Calculate the decision boundary line
    plot_y = (-1./theta(3)).*(theta(2).*plot_x + theta(1));


    plot(plot_x, plot_y);

% ------------------- 2 b over ----------------------

%--------------------- 4a --------------------------- 
X = load('q4x.dat');
tempy = importdata('q4y.dat');

m = length(tempy);
y = zeros(m,1);
for i=1:m
	if(strcmp(tempy(i),'Canada'))
		y(i) = 1;
	
end
end

%adding the intercept term
%X = [ones(m,1) X];

n = size(X,2);
pos = find(y==1); neg = find(y==0);
% --------------- calculating phi -------------------

phi = sum(y(pos,:));

phi = phi/m;

% --------------- calculated phi --------------------

% --------------- calculating mu0 -------------------

mu0 = zeros(n,1);

mu0 = transpose(sum(X(neg,:))./length(y(neg,:)));


% --------------- calculated mu0 -------------------

% --------------- calculating mu1 -------------------

mu1 = zeros(n,1);

mu1 = transpose(sum(X(pos,:))./length(y(pos,:)));


% --------------- calculated mu1 -------------------

% --------------- calculating covariance matrix sigma -------------------

covmat = zeros(n,n);

for i=1:m
	if(y(i,:)==0)
		mu=mu0;
	else
		mu=mu1;
	end
	covmat = covmat + (transpose(X(i,:))-mu)*transpose(transpose(X(i,:))-mu);
	end;

	covmat = covmat./m;

% --------------- calculated covariance matrix sigma -------------------

plot(X(pos,1),X(pos,2),'+');
hold
plot(X(neg,1),X(neg,2),'o');

plo_x = [min(X(:,1))-2,  max(X(:,1))+2];

c1 = pinv(covmat)*(mu1 - mu0);
c2 = (transpose(mu1) - transpose(mu0))*pinv(covmat);
c3 = transpose(mu0)*pinv(covmat)*mu0;
c4 = transpose(mu1)*pinv(covmat)*mu1;

plo_y = (c4 - c3 - (c2(1)*plo_x) - (c1(1)*plo_x))/(c1(2)+c2(2)) + log(exp(0.5))

plot(plo_x, plo_y);

% --------------------- 4th second part ----------------------

covmat0 = zeros(n,n);

for i=1:m
	if(y(i,:) == 0)
		covmat0 = covmat0 + (transpose(X(i,:))-mu0)*transpose(transpose(X(i,:))-mu0);

end;
end;
covmat0 = covmat0 ./ m;


covmat1 = zeros(n,n);
for i=1:m
	if(y(i,:) == 1)
		covmat1 = covmat1 + (transpose(X(i,:))-mu1)*transpose(transpose(X(i,:))-mu1);

end;
end;
covmat1 = covmat1 ./ m;

c2 = pinv(covmat1)*mu1 - pinv(covmat0)*mu0;
c3 = transpose(mu1)*pinv(covmat1) - transpose(mu0)*pinv(covmat0);
c4 = (transpose(mu0)*pinv(covmat0)*mu0) - (transpose(mu1)*pinv(covmat1)*mu1);



syms x0;
syms x1;
x = [x0;x1];

soln = solve((((transpose(x-mu0))*(pinv(covmat0))*(x-mu0)) - ((transpose(x-mu1))*(pinv(covmat1))*(x-mu1))),x1);
soln = soln(1);
x0 = linspace(40,200,100);
plot(x0,eval(soln));








