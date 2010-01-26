function engine = global_joint_inf_engine(bnet)
% GLOBAL_JOINT_INF_ENGINE Construct the global joint distribution as a potential
% engine = global_joint_inf_engine(bnet)
%
% Warning: this has size exponential in the number of discrete hidden variables

engine.jpot = [];
engine = class(engine, 'global_joint_inf_engine', inf_engine(bnet));    
