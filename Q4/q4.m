X = load('q4x.dat');
tempy = importdata('q4y.dat');

m = length(tempy);
y = zeros(m,1);
for i=1:m
	if(strcmp(tempy(i),'Canada'))
		y(i) = 1;
	
end
end

n = size(X,2);
pos = find(y==1); neg = find(y==0);

%---------------------- phi ---------------------
phi = sum(y(pos,:));
phi = phi/m;
%---------------------- phi ---------------------


%------------------------------------ a --------------------------------------------
mu0 = zeros(n,1);
mu0 = transpose(sum(X(neg,:))./length(y(neg,:)));

mu1 = zeros(n,1);
mu1 = transpose(sum(X(pos,:))./length(y(pos,:)));


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
%-----------------------------------------------------------------------------------

disp('values of mu0, mu1 and covariance matrix are as follows....');
mu0
mu1
covmat
disp('press enter to continue..');
pause;




%--------------------- b -----------------------------
disp('plotting graphs');
pause(1);
hold off;
plot(X(pos,1),X(pos,2),'+');
hold
plot(X(neg,1),X(neg,2),'o');
xlabel('growth ring diameter in fresh water');
ylabel('growth ring diameter in marine water');
legend('Canada', 'Alaska');
%-----------------------------------------------------


%--------------------------------------------- c ------------------------------------------
plot_x = [min(X(:,1))-0.1,  max(X(:,1))+0.1];

c1 = pinv(covmat)*(mu1 - mu0);
c2 = (transpose(mu1) - transpose(mu0))*pinv(covmat);
c3 = transpose(mu0)*pinv(covmat)*mu0;
c4 = transpose(mu1)*pinv(covmat)*mu1;

plot_y = (c4 - c3 - (c2(1)*plot_x) - (c1(1)*plot_x) + (2*log((1-phi)/phi)) )/(c1(2)+c2(2));

disp(' press enter to plot the boundary line ');
pause;
plot(plot_x, plot_y);
%------------------------------------------------------------------------------------------
hold off;
disp('press enter to continue to part d .....');
pause;

%------------------------------------------ d -----------------------------------------------
covmat0 = zeros(n,n);

for i=1:m
	if(y(i,:) == 0)
		covmat0 = covmat0 + (transpose(X(i,:))-mu0)*transpose(transpose(X(i,:))-mu0);
    end;
end;
covmat0 = covmat0 ./ length(y(neg,:));


covmat1 = zeros(n,n);
for i=1:m
	if(y(i,:) == 1)
		covmat1 = covmat1 + (transpose(X(i,:))-mu1)*transpose(transpose(X(i,:))-mu1);
    end;
end;
covmat1 = covmat1 ./ length(y(pos,:));

disp('the values of mu0, mu1, covariance matrix 0 and covariance matrix 1 are as follows');
mu0
mu1
covmat0
covmat1
%--------------------------------------------------------------------------------------------

%---------------------------------------------- e -------------------------------------------
syms x0;
syms x1;
x = [x0;x1];

%soln = solve(((log((1-phi)/phi))*(sqrt(det(covmat1))/sqrt(det(covmat0))))+(-0.5)*(((transpose(x-mu0))*(pinv(covmat0))*(x-mu0)) - ((transpose(x-mu1))*(pinv(covmat1))*(x-mu1))),x1);
soln = solve((log(((1-phi)/phi)*((sqrt(det(covmat1)))/(sqrt(det(covmat0))))))+(-0.5)*(((transpose(x-mu0))*(pinv(covmat0))*(x-mu0)) - ((transpose(x-mu1))*(pinv(covmat1))*(x-mu1))),x1);




%soln = soln(2);
x0 = linspace(min(X(:,1))-0.1,  max(X(:,1))+0.1,100);
temp1 = eval(soln(1));
if(temp1(1)<0)
soln = soln(2);
else
soln = soln(1);
end;

disp('plotting graphs ...');
pause(1);
figure();
plot(X(pos,1),X(pos,2),'+');
hold;
plot(X(neg,1),X(neg,2),'o');
xlabel('growth ring diameter in fresh water');
ylabel('growth ring diameter in marine water');
legend('Canada', 'Alaska');
plot(x0,eval(soln));
%--------------------------------------------------------------------------------------------















