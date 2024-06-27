function [h_approx] = policy_fun(theta,X,model,policy)
  x_dim = model.x_dim;
  u_dim = 1;
  neurons = policy.approx.neurons;
  W = reshape(theta(1:(x_dim+1)*neurons),(x_dim+1),neurons);
  C = reshape(theta((x_dim+1)*neurons+1:end),neurons,u_dim);
  h_aprox = mlp_out(C,W,X,model,policy);
end
