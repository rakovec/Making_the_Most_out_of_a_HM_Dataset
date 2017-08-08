function [theta_hathat,theta_BS]=strapest(nboot,bootfun,varargin)
% STRAPEST biased-corrected Bootstrap estimator.

theta_hat=feval(bootfun,varargin{:});
theta_BS=feval(@bootstrp,nboot,bootfun,varargin{:});
theta_hathat=2*theta_hat-theta_BS;
end
