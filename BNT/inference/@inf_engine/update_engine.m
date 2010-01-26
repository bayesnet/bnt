function engine = update_engine(engine, newCPDs)
% UPDATE_ENGINE Update the engine to take into account the new parameters (inf_engine).
% engine = update_engine(engine, newCPDs)
%
% This generic method is suitable for engines that do not process the parameters until 'enter_evidence'.

engine.bnet.CPD = newCPDs;
