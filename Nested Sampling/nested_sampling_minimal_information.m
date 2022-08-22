%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimal information model implementation of the Ornstein-Uhlenbeck process
% with different characteristic time values per time interval n with the 
% amount of time intervals also being a parameter.
%
% Contributors to the programming: Michael Lomholt, Maximilian Konrad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs=importdata('Patient_I.csv'); % Import the data

% Define inverse prior
n_prior=@(out_wish) ns_exp_disc(Inf,5,10,out_wish); % Prior for time interval is discrete uniform between 5 and 10
n_prob=n_prior('p');
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
sigma_invprior=@(u) 0.01*(10.^(1*inv_normal(u)));
tau_invprior=@(u) 50*(10.^(1*inv_normal(u)));

% Define model
model=struct;
model.genu=@(obs) genu_func(n_prior('n'));
model.adjust_u=@(u,obs) adjust_u_func(u,n_prior('n'));
model.invprior=@(u,obs) invprior_func(u,obs,n_prior('n'),tau_invprior,sigma_invprior);
model.logl=@(obs,theta) log_likelihood_multi_tau(length(obs),theta.x_hat,theta.t_switch,theta.tau,theta.sigma_x,obs);
model.prior_disc=@(disc) n_prob(disc);
model.cont=@(theta) [theta.x_hat theta.t_switch theta.tau theta.sigma_x];
model.disc=@(theta) [theta.n]
model.labels=@(disc,obs) label_func(disc);
model.names=@(disc) sprintf('%i intervals',disc);

% Define algorithm
options=struct;
options.algorithm=@ns_algorithm;
misc=struct;

% Save results
results=cube_main(obs,model,misc,options);

% Generate u value
% Necessary for accounting for changing n
function u = genu_func(n_invprior)
  u.n=rand;
  n=n_invprior(u.n);
  u.x_hat=rand(1,n);
  u.t_switch=rand(1,n-1);
  u.tau=rand(1,n);
  u.sigma_x=rand;
end

% Adjust u function
% Necessary for accounting for changing n
function uout = adjust_u_func(u,n_invprior,dim)
  if isfield(u,'n') && length(u.n)>0
    uout.n=u.n(1);
  else
    uout.n=rand;
  end
  n=n_invprior(uout.n);
  fields={'x_hat','t_switch','tau','sigma_x'};
  expected=[n n-1 n 1];
  uout=ns_adjust(u,uout,fields,expected);
end

% Define thetas in function
% Necessary for accounting for changing n
function theta = invprior_func(u,obs,n_invprior,tau_invprior,sigma_invprior)
  theta.n=n_invprior(u.n);
  theta.x_hat=cumsum(pf_inv_simple_dirichlet(u.x_hat,4)*(max(obs)-min(obs)))+min(obs);
  theta.t_switch=cumsum(pf_inv_simple_dirichlet(u.t_switch,4)*(length(obs)-1-theta.n)+1);
  theta.tau=tau_invprior(u.tau);
  theta.sigma_x=sigma_invprior(u.sigma_x);
end

% Label results in result file
function label = label_func(disc)
  n=disc;
  label1={};
  label2={};
  label3={};
  for j=1:n
    label1=[label1 {sprintf('#%i x_hat: ',j)}];
    label3=[label3 {sprintf('#%i tau: ', j)}];
    if j~=n
      label2=[label2 {sprintf('#%i t_switch: ',j)}];
    end
  end
  label=[label1 label2 label3 {'sigma_x:'}];
end

