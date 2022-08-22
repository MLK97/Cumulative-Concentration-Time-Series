%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-tau approximation implementation of the Ornstein-Uhlenbeck process
% with the same characteristic time value over all time intervals.
%
% Contributors to the programming: Michael Lomholt, Maximilian Konrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%obs=load('data.txt');
%data=importdata('data.csv');
data=importdata('data_five_intervals.csv');
obs=data.data(1:1000);

n=5;
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
sigma_invprior=@(u) 0.01*(10.^(1*inv_normal(u)));
tau_invprior=@(u) 50*(10.^(1*inv_normal(u)));
model=struct;
model.genu=@(obs) struct('x_hat',rand(1,n),'t_switch',rand(1,n-1),'tau',rand,'sigma_x',rand);
model.invprior=@(u,obs) struct('x_hat',cumsum(pf_inv_simple_dirichlet(u.x_hat,4)*(max(obs)-min(obs)))+min(obs),...
                               't_switch',cumsum(pf_inv_simple_dirichlet(u.t_switch,4)*(length(obs)-1-n)+1),...
                               'tau',tau_invprior(u.tau),...
                               'sigma_x',sigma_invprior(u.sigma_x));
model.logl=@(obs,theta) log_likelihood(length(obs),theta.x_hat,theta.t_switch,theta.tau,theta.sigma_x,obs);
model.cont=@(theta) [theta.x_hat theta.t_switch theta.tau theta.sigma_x];
model.labels=@(disc,obs) label_func(disc,n);

options=struct;
options.algorithm=@ns_algorithm;

results=cube_main(obs,model,struct,options);

function label = label_func(disc,n)
  label1={};
  label2={};
  for j=1:n
    label1=[label1 {sprintf('#%i x_hat: ',j)}];
    if j~=n
      label2=[label2 {sprintf('#%i t_switch: ',j)}];
    end
  end
  label=[label1 label2 {'tau:'} {'sigma_x:'}];
end

