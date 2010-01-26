function engine = set_fwdback(engine, fb)
% SET_FWDBACK Set the field 'fwdback', which contains the frontiers after propagation
% engine = set_fwdback(engine, fb)
%
% This is used by frontier_fast_inf_engine/enter_evidence
% as a workaround for Matlab's annoying privacy control
    
engine.fwdback = fb;
