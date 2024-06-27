function out = mlp_out(C,W,x,model,policy)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% x_num = size(x,2);
% out = C'*[tanh(W'*[x;ones(1,x_num)]);ones(1,x_num)];


switch policy.approx.type
    case 'NN1'
        x_num = size(x,2);
        out = C'*[tanh(W'*[x;ones(1,x_num)]);ones(1,x_num)];
    case 'NN2'
        x_num = size(x,2);
        out = C'*tanh(W'*[x;ones(1,x_num)]);
    case 'NN3'
        out = C'*tanh(W'*x);
end

% out = model.u1_bounds(2)*tanh(0.5*out); %limited output
% out = model.u1_bounds(2)*tanh(0.5*out); %limited output


end

