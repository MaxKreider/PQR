%we consider the model
%
%  dx = [omega + epsilon*cos(x) + alpha*sin(t+phi) + beta*sin(2*t)]dt
%           + sqrt(2D)dW(t)
%
%and use the matrix continued fraction method to solve for the
%cyclostationary distribution.

%% setup

%mr clean
clc

%cutoff point to approximate infinite continued fraction as a finite one
cutoff = 20;
cutoff = cutoff*2+1;

%mesh for eigenfunction
dt = 0.05;
X = 0:dt:2*pi;
Y = X';

%system parameters
D = .01;
omega = 1;
epsilon = .2;
alpha = .5;
beta = .5;

%cycle through phase offset
phi = 0:0.1:2*pi;

%firing rate solution matrix
J_matrix = zeros(length(phi), length(Y));


%% explicit matrix method

%loop over phase offsets
for pp=1:length(phi)

    %initialize matrix
    my_matrix_friend = zeros(cutoff^2,cutoff^2);

    %insert small matrices into big matrix
    count = 0;
    jcount = -(cutoff-1)/2;
    for j=0:cutoff-1

        %block diagonal
        my_matrix_friend(count+1:cutoff*(j+1),count+1:cutoff*(j+1)) = Q0(jcount,D,omega,epsilon,alpha,beta,phi(pp),cutoff);

        %block subdiagonal
        if (j >= 1)
            my_matrix_friend(count+1:cutoff*(j+1),count-cutoff+1:cutoff*(j)) = Qm(jcount,D,omega,epsilon,alpha,beta,phi(pp),cutoff);
        end

        %block superdiagonal
        if (j < cutoff-1)
            my_matrix_friend(count+1:cutoff*(j+1),count+cutoff+1:cutoff*(j+2)) = Qp(jcount,D,omega,epsilon,alpha,beta,phi(pp),cutoff);
        end

        %update count
        count = count + cutoff;
        jcount = jcount + 1;
    end

    %diagonalize
    [xi,lambda_boi] = eig(my_matrix_friend);
    lambda_boi = diag(lambda_boi);

    %sort
    [~,temp] = sort(abs(lambda_boi));
    lambda_boi = lambda_boi(temp);
    xi = xi(:,temp);

    %display first 10 eigenvalues of smallest magnitude
    lambda_boi(1:30);


    %% construct cyclostationary distribution in terms of basis functions

    %select eigenvector of choice
    index = 1;

    %select eigenvector of choice
    xi = xi(:,index)/max(xi(:,index));

    %evaluate
    cyclo = zeros(length(X),length(Y));
    basis = zeros(cutoff,cutoff);
    for p=1:length(X)
        for q=1:length(Y)

            %form the basis at each point in the domain
            for j=-(cutoff-1)/2:(cutoff-1)/2
                for k=-(cutoff-1)/2:(cutoff-1)/2
                    basis(j+(cutoff-1)/2+1,k+(cutoff-1)/2+1) = exp(sqrt(-1)*(j*X(p)+k*Y(q)));
                end
            end
            new_basis = reshape(basis',[cutoff^2, 1]);

            %insertion
            cyclo(p,q) = dot(xi,new_basis);
        end
    end

    %Matlab indexing is dumb, and normalize
    cyclo = cyclo.';
    cyclo = cyclo/(sum(sum(cyclo)*dt*dt));

    %prune imaginary part
    cyclo = real(cyclo);


    %% visualize eigenvalues (pruned)

    % %plot the eigenvalues
    % figure(1)
    % plot(lambda_boi(1:31),'k.','markersize',40)
    % box on
    % grid on
    % xlabel('Re(\lambda)')
    % ylabel('Im(\lambda)')
    % set(gca,'fontsize',15)
    %
    % %plot it
    % figure(2)
    % imagesc(X,Y,cyclo)
    % box on
    % axis square
    % set(gca,'fontsize',15)
    % xlabel('x')
    % ylabel('t')
    % colorbar
    % colormap hot
    % title('Cyclostationary distribution')
    % set(gca,'Ydir','normal')


    %% extract firing rate, i.e., J(x,t) at x = 2*pi

    %compute x-derivative
    [dx,~] = gradient(cyclo, dt);

    %evaluate drift term at x = 2pi
    x_val = 2*pi;
    a_t = omega + epsilon*cos(x_val) + alpha*sin(Y + phi(pp)) + beta*sin(2*Y);

    %cyclo at the edge
    cyclo_edge = cyclo(:, end);
    der_edge = dx(:, end);

    %flux
    J = a_t .* cyclo_edge - D * der_edge;
    
    figure(4)
    hold on
    plot3(Y, phi(pp)*ones(size(Y)), J, 'k-', 'LineWidth', 1.5)
    xlabel('time t')
    ylabel('phase offset \phi')
    zlabel('firing rate J(t;\phi)')
    title('Firing rate')
    view(45, 30)
    box on
    set(gca, 'fontsize', 12)

    figure(5)
    hold on
    plot(phi(pp),mean(J),'k.','markersize',40)
    xlabel('phase offset \phi')
    ylabel('E(J)')
    title('Mean firing rate')
    box on
    axis square
    set(gca,'fontsize',12)

    figure(6)
    hold on
    plot(phi(pp),max(J),'k.','markersize',40)
    xlabel('phase offset \phi')
    ylabel('max(J)')
    title('Max firing rate')
    box on
    axis square
    set(gca,'fontsize',12)

    pp

end


%% functionals to make my life easier

%Q-
function[u] = Qm(j,D,omega,epsilon,alpha,beta,phi,cutoff)

%initialize matrix
u = zeros(cutoff,cutoff);

%define the k interval
interval = (cutoff-1)/2;

%fill him up
count = 1;
for k=-interval:interval
    u(count,count) = -epsilon/2*sqrt(-1)*j;
    count = count+1;
end

%make sure that the size is correct
u = u(1:cutoff,1:cutoff);
end

%Q
function[u] = Q0(j,D,omega,epsilon,alpha,beta,phi,cutoff)

%initialize matrix
u = zeros(cutoff,cutoff);

%define the k interval
interval = (cutoff-1)/2;

%fill him up
count = 1;
for k=-interval:interval
    u(count,count+2) = beta*(j/2);
    u(count,count+1) = alpha/2*exp(-sqrt(-1)*phi)*(j);
    u(count,count) = -D*j^2-sqrt(-1)*(k+omega*j);
    u(count+1,count) = -alpha/2*exp(sqrt(-1)*phi)*(j);
    u(count+2,count) = -beta*(j/2);
    count = count+1;
end

%make sure that the size is correct
u = u(1:cutoff,1:cutoff);
end

%Q+
function[u] = Qp(j,D,omega,epsilon,alpha,beta,phi,cutoff)

%initialize matrix
u = zeros(cutoff,cutoff);

%define the k interval
interval = (cutoff-1)/2;

%fill him up
count = 1;
for k=-interval:interval
    u(count,count) = -epsilon/2*sqrt(-1)*j;
    count = count+1;
end

%make sure that the size is correct
u = u(1:cutoff,1:cutoff);
end





