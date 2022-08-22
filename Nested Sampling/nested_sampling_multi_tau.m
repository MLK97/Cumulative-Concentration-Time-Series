%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal Model implementation of the Ornstein-Uhlenbeck process
% with a individual characteristic time per time interval.
%
% Contributors to the programming: Michael Lomholt, Maximilian Konrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs=importdata('Patient_A.csv').data; % Import the data
n=9; % Define the number of time intervals

% Define inverse priors
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
sigma_invprior=@(u) 0.01*(10.^(1*inv_normal(u)));
tau_invprior=@(u) 150*(10.^(1*inv_normal(u)));

% Define model
model=struct;
model.genu=@(obs) struct('x_hat',rand(1,n),'t_switch',rand(1,n-1),'tau',rand(1,n),'sigma_x',rand);
model.invprior=@(u,obs) struct('x_hat',cumsum(pf_inv_simple_dirichlet(u.x_hat,4)*(max(obs)-min(obs)))+min(obs),...
                               't_switch',cumsum(pf_inv_simple_dirichlet(u.t_switch,4)*(length(obs)-1-n)+1),...
                               'tau',tau_invprior(u.tau),...
                               'sigma_x',sigma_invprior(u.sigma_x));

model.logl=@(obs,theta) log_likelihood_multi_tau(length(obs),theta.x_hat,theta.t_switch,theta.tau,theta.sigma_x,obs);
model.cont=@(theta) [theta.x_hat theta.t_switch theta.tau theta.sigma_x];
model.labels=@(disc,obs) label_func(disc,n);

% Define algorithm
options=struct;
options.algorithm=@ns_algorithm;

% Save results
results=cube_main(obs,model,struct,options);

% Label results in result file
function label = label_func(disc,n)
  label1={};
  label2={};
  label3={};
  for j=1:n
    label1=[label1 {sprintf('#%i x_hat: ',j)}];
    label3=[label3 {sprintf('#%i tau: ',j)}];
    if j~=n
      label2=[label2 {sprintf('#%i t_switch: ',j)}];
    end  
  end
  label=[label1 label2 label3 {'sigma_x:'}];
end
