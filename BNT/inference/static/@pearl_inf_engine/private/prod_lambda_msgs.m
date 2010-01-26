function lam = prod_lambda_msgs(n, cs, msg, msg_type, except)

if nargin < 5, except = -1; end

lam = msg{n}.lambda_from_self;
switch msg_type
  case 'd',
   for i=1:length(cs)
     c = cs(i);
     if c ~= except
       lam = lam .* msg{n}.lambda_from_child{i};
     end
   end  
 case 'g',
  if isinf(lam.precision) % isfield(lam, 'observed_val')
    return; % pass on the observed msg
  end
   for i=1:length(cs)
     c = cs(i);
     if c ~= except
       m = msg{n}.lambda_from_child{i};
       lam.precision = lam.precision + m.precision;
       lam.info_state = lam.info_state + m.info_state;
     end
   end  
end



