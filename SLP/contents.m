% BNT Structure Learning Package 
% http://banquiseasi.insa-rouen.fr/projects/bnt-slp/
%
% Some usefull add-ons for Bayes Net Toolbox
% ==================================================================================================
% v1.5 : 20 avril 2008
% ==================================================================================================
%
% Data manipulation
% -----------------
%
%  - mat_to_bnt             :   matrix to cell array conversion (with missing data encoding in the matrix)
%  - hist_ic                :   optimal Histogram based on IC information criterion (PhL)
%  - histc_ic               :   Histogram count , for HIST_IC edges (PhL)
%
%
% Bayes Net Structure Learning
% ----------------------------
%
%  - cond_indep_chisquare   :   test if X indep Y given Z using ChiSquare test (Pearson's or Likelihood Ratio Test) (PhL,OF)
%                               (for use with LEARN_STRUCT_PDAG_PC or LEARN_STRUCT_PDAG_IC_STAR)
%  - test_chisquare         :   test for COND_INDEP_CHISQUARE function
%  - test_pc                :   test for LEARN_STRUCT_PDAG_PC function with ChiSquare test (PhL)
%
%  - mk_alarm_bnet	    :   make the bnet of ALARM network (WH) (used in test_sem3)
%  - mk_asia_bnet	    :   make the bnet of ASIA network (PhL,OF) (used in test functions)
%
%  - mk_naive_struct        :   generate the naive bayes structure (for a given class node)
%
%  - learn_struct_mwst      :   structure learning using maximum spanning tree (OF, PhL)
%  - mutual_info_score      :   mutual information scoring (OF, PhL, WXY)
%  - test_mwst              :   test for LEARN_STRUCT_MWST function (PhL)
%
%  - learn_struct_tan       :   structure learning giving the best tree augmented naive bayes classifier (OF, PHL, NS)
%
%  - learn_struct_hc        :   structure learning using hill climbing (GL)
%  - learn_struct_gs        :   structure learning using greedy search (GL)
%  - learn_struct_gs2       :   structure learning using greedy search with cache implementation (GL, WH, OF, PhL)
%                               (use mk_nbrs_of_dag_topo (developped WH) instead of mk_nbrs_of_dag)
%  - score_init_cache       :   cache initialisation for local score computation (OF, PhL) 
%  - score_family           :   new version of BNT function with cache implementation (OF, PhL)
%  - score_dags             :   new version of BNT function with cache implementation (OF, PhL, DH)
%  - test_gs2               :   test for LEARN_STRUCT_GS2 function (PhL)
%
%  - learn_struct_ges	    :   structure learning using Greedy Equivalence Search (PhL)
%  - mk_nbrs_of_pdag_add    :   generate the sup. inclusion boundary of a given pdag (PhL)
%  - mk_nbrs_of_pdag_del    :   generate the inf. inclusion boundary of a given pdag (PhL)
%  - test_ges               :   test for LEARN_STRUCT_GES function (PhL)
%
%  - learn_struct_EM        :   structure learning using structural EM (WH)
%  - multiply_one_marginal.c
%  - mk_nbrs_of_dag_topo    :   generate the neighbours of a given dag (WH)
%                               (better implementation than mk_nbrs_of_dag)
%  - test_sem1              :   demo 1 (SPRINKER)
%  - test_sem2              :   demo 2 (DISCRETE1)
%  - test_sem3              :   demo 3 (ALARM)
%
%  - learn_struct_mwst_EM   :   MWST structure learning with missing data (OF, PhL)
%
%  - cpdag_to_dag           :   return a dag for a given CPDAG (OF, PhL)
%  - dag_to_cpdag           :   return the CPDAG, representant of the equivalent class of the dag (OF, PhL)
%  - pdag_to_dag            :   return a DAG that instantiates the given pdag [Dor&Tarsi] (PhL, OF)
%  - test_cpdag             :   test CPDAG_to_DAG and DAG_to_CPDAG functions (PhL)
%  - kl_divergence          :   Kullback-Leibler divergence between two bnet distributions (PhL)
%  - kl_divergence2         :   Kullback-Leibler divergence between two bnet distributions (PhL) 
%
%  - test_structure         :   runs all the test functions (PhL)
%
%  - gener_MCAR_net         :   to gener a BN which modelise a process of incomplete data generation with MCAR assumptions (OF)
%  - gener_MAR_net          :   to gener a BN which modelise a process of incomplete data generation with MAR assumptions (OF)
%  - gener_data_from_bnet_miss: to gener a incomplete dataset from a BN creted with gener_MCAR_net or gener_MAR_net (OF)
%
%  - editing_dist           :   editing distance between two DAG (PhL)
%                               (memory optimization, but quite slow !)
%  - learn_struct_bnpc      :   structure learning with BN-Power constructor (OF, PhL)
%  - learn_struct_tan_EM    :   Tree Augmented Naive Bayes structure learning with missing data (OF, PhL)
%  - learn_struct_ges_EM    :   structure learning in Markov equivalent space with missing data (OF, HB, PhL)
%  - ......
%
%
% Contributors :
% ------------
% LITIS Rouen, France
%   - PhL   : Philippe Leray (philippe.leray@univ-nantes.fr)
%   - OF    : Olivier Francois (francois.olivier.c.h@gmail.com)
%
% External parts
%   - GL    : Gang Li, Deakin University (gangli@deakin.edu.au)
%   - WH    : Wei Hu, Intel (wei.hu@intel.com)
%   - DH    : Derek Hoiem (dhoiem@cs.cmu.edu)
%   - NS    : Navid Serrano, Jet Propulsion Laboratory (Navid.Serrano@jpl.nasa.gov)
%   - WXY   : Wang Xiang Yang, Shanghai JiaoTong University (wangxiangyang@sjtu.edu.cn)
%   - HB    : Hanene Borchani (hanene.borchani@gmail.com)
%

