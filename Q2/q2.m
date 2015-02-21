hold off;
X = load('q2x.dat');
y = load('q2y.dat');

m = length(y);
%adding the intercept term
X = [ones(m,1) X];

theta = zeros(size(X,2),1);

n = length(theta);

disp(' working out part a....');
pause(1);

%------------------------- a ---------------------------
oldtheta = 100.*ones(n,1);
epsilon = 0.01;

while(norm(theta - oldtheta)>epsilon)
	oldtheta = theta;
	theta = theta - pinv(hessian(X,y,theta))*grad(X,y,theta);
end;

disp(' resulting cofficients are..');
theta
disp(' press enter to continue');
pause; 

%-------------------------------------------------------

disp(' generating graphs and boundary line.... ');
pause(1);

%------------------------- b ----------------------------

pos = find(y==1); neg = find(y==0);

plot(X(pos,2),X(pos,3),'+');
hold
plot(X(neg,2),X(neg,3),'o');

xlabel('x1');
ylabel('x2');

% extracting two points to define the line 
    plot_x = [min(X(:,2))-0.1,  max(X(:,2))+0.1];

% Calculate the decision boundary line
    plot_y = (-1./theta(3)).*(theta(2).*plot_x + theta(1));


    plot(plot_x, plot_y);
    legend('y=1', 'y=0');

%-------------------------------------------------------

hold off;
disp(' press enter to continue ');
pause;



