function [engine, loglik] = enter_evidence(engine, evidence, varargin)
% ENTER_EVIDENCE enter evidence to engine including discrete and continuous evidence
% [engine, ll] = enter_evidence(engine, evidence)
%
% ll is always 0, which is wrong.

if ~isempty(engine.evidence)
    bnet = bnet_from_engine(engine);
    engine = stab_cond_gauss_inf_engine(bnet);
    engine.evidence = evidence;
else
    engine.evidence = evidence;
    bnet = bnet_from_engine(engine);
end

engine.evidence = evidence;
bnet = bnet_from_engine(engine);
ns = bnet.node_sizes(:);
observed = ~isemptycell(evidence);
onodes = find(observed);
hnodes = find(isemptycell(evidence));
cobs = myintersect(bnet.cnodes, onodes);
dobs = myintersect(bnet.dnodes, onodes);

engine = incorporate_dis_evidence(engine, dobs, evidence);
l = length(cobs);
for i = 1:l
    node = cobs(i);
    engine = incorporate_singleconts_evidence(engine, node, evidence);
end
clpot = engine.clpot;

clq_num = length(engine.cliques);
for n=engine.postorder(1:end-1)
  for p=parents(engine.jtree, n)
      [margpot, comppot] = complement_pot(clpot{n}, engine.separator{p,n});
      clpot{n} = comppot;
      clpot{p} = combine_pots(clpot{p}, margpot);
  end
end

temppot = clpot;
for n=engine.preorder
  for c=children(engine.jtree, n)
    seppot{n,c} = marginalize_pot(temppot{n}, engine.separator{n,c});
    temppot{c} = direct_combine_pots(temppot{c}, seppot{n,c});
  end
end
engine.clpot = clpot;
engine.seppot = seppot;

[pot,loglik]=normalize_pot(clpot{engine.root});

%%%%%%%%%%%%%%%%%%
function engine = incorporate_dis_evidence(engine, donodes, evidence)
l = length(donodes);
for i=donodes(:)'
    node = i;
    clqid = engine.clq_ass_to_node(node);
    pot = struct(engine.clpot{clqid});
    ns = zeros(1, max(pot.domain));
    ns(pot.ddom) = pot.dsizes;
    ns(pot.cheaddom) = pot.cheadsizes;
    ns(pot.ctaildom) = pot.ctailsizes;
    ddom = pot.ddom;
    
    potcarray = cell(1, pot.dsize);
    for j =1:pot.dsize
        tpotc = struct(pot.scgpotc{j});
        potcarray{j} = scgcpot(tpotc.cheadsize, tpotc.ctailsize, 0, tpotc.A, tpotc.B, tpotc.C);
    end
    
    if length(ns(ddom)) == 1
        matrix = pot.scgpotc;
    else
        matrix = reshape(pot.scgpotc,ns(ddom)); 
        potcarray = reshape(potcarray, ns(ddom));
    end
    
    map = find_equiv_posns(node, ddom);
    vals = cat(1, evidence{node});
    index = mk_multi_index(length(ddom), map, vals);
    potcarray(index{:}) = matrix(index{:});
    potcarray = potcarray(:);
    %keyboard;
    engine.clpot{clqid} = scgpot(pot.ddom, pot.cheaddom, pot.ctaildom, ns, potcarray);
end

%%%%%%%%%%%%%%%%%%
function engine = incorporate_singleconts_evidence(engine, node, evidence)
%incorporate_singleconts_evidence incorporate evidence of 1 continuous node
B = engine.cliques_bitv;
clqs_containnode = find(all(B(:,node), 2)); % all selected columns must be 1
% Every continuous node necessarily apears as head in exactly one clique,
% which is the clique where it appears closest to the strong root. In all other
% clique potentials where it appears, it must be a tail node.
clq_ev_as_head = [];
for i = clqs_containnode(:)'
    pot = struct(engine.clpot{i});
    if myismember(node, pot.cheaddom)
        clq_ev_as_head = [clq_ev_as_head i];
        break;
    end
end
	       
% If we will incorporate the evidence node which is head of a potential we must rearrange
% the juntion tree by push operation until the tail of the include potential is empty
if ~isempty(clq_ev_as_head)
    assert(1 == length(clq_ev_as_head));
    i = clq_ev_as_head;
    pot = struct(engine.clpot{i});
    while ~isempty(pot.ctaildom)
        [engine, clqtoroot] = push(engine, i, node);
        i = clqtoroot;
        pot = struct(engine.clpot{i});
    end
    B = engine.cliques_bitv;
    clqs_containnode = find(all(B(:,node), 2));
end

for i = clqs_containnode(:)'
    pot = struct(engine.clpot{i});
    if myismember(node, pot.cheaddom)
        engine.clpot{i} = incoporate_evidence_headnode(engine.clpot{i}, node, evidence);
    else
        %assert(myismember(node, pot.ctaildom));
        engine.clpot{i} = incoporate_evidence_tailnode(engine.clpot{i}, node, evidence);
    end
end

%%%%%%%%%%%%%%%%%%
function newscgpot = incoporate_evidence_tailnode(pot, node, evidence)
%ENTER_EVIDENCE_TAILNODE enter the evidence of 1 tailnode of the scgpot
newscgpot = pot;
pot = struct(pot);
%if isempty(pot.ctaildom)
if ~myismember(node, pot.ctaildom)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In this case there is no real dependency of the head nodes %
    % on the tail. The potential should be returned unchanged    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    return;
end
%newscgpot = scgpot([], [], [], []);
assert(myismember(node, pot.ctaildom));
ni = block(find_equiv_posns(node, pot.ctaildom), pot.ctailsizes);

ctaildom = mysetdiff(pot.ctaildom, node);
cheaddom = pot.cheaddom;
ddom = pot.ddom;
domain = mysetdiff(pot.domain, node);
dsize = pot.dsize;
ns = zeros(1, max(pot.domain));
ns(pot.ddom) = pot.dsizes;
ns(pot.cheaddom) = pot.cheadsizes;
ns(pot.ctaildom) = pot.ctailsizes;
cheadsizes = pot.cheadsizes;
cheadsize = pot.cheadsize;
ctailsizes = ns(ctaildom);
ctailsize = sum(ns(ctaildom));

potarray = cell(1, dsize);
for i=1:dsize
    potc = struct(pot.scgpotc{i});
    B = potc.B;
    A = potc.A + B(:, ni)*evidence{node};
    B(:, ni) = [];
    potarray{i} = scgcpot(cheadsize, ctailsize, potc.p, A, B, potc.C);
end

newscgpot = scgpot(ddom, cheaddom, ctaildom, ns, potarray);

%%%%%%%%%%%%%%%%
function newscgpot = incoporate_evidence_headnode(pot, node, evidence)
%ENTER_EVIDENCE_HEADNODE 
pot = struct(pot);
y2 = evidence{node};
assert(myismember(node, pot.cheaddom));
assert(isempty(pot.ctaildom));
ddom = pot.ddom;
cheaddom = mysetdiff(pot.cheaddom, node);
ctaildom = pot.ctaildom;
dsize = pot.dsize;
domain = mysetdiff(pot.domain, node);

ns = zeros(1, max(pot.domain));
ns(pot.ddom) = pot.dsizes;
ns(pot.cheaddom) = pot.cheadsizes;
ns(pot.ctaildom) = pot.ctailsizes;
ctailsizes = ns(ctaildom);
ctailsize = sum(ctailsizes);
cheadsizes = ns(cheaddom);
cheadsize = sum(cheadsizes);
onodesize = ns(node);

p = zeros(1,dsize);
A1 = zeros(cheadsize, dsize);
A2 = zeros(onodesize, dsize);
C11 = zeros(cheadsize, cheadsize, dsize);
C12 = zeros(cheadsize, onodesize, dsize);
C21 = zeros(onodesize, cheadsize, dsize);
C22 = zeros(onodesize, onodesize, dsize);
ZM = zeros(onodesize, onodesize);

n1i = block(find_equiv_posns(cheaddom, pot.cheaddom), pot.cheadsizes);
n2i = block(find_equiv_posns(node, pot.cheaddom), pot.cheadsizes);

indic = 0;
for i=1:dsize
    potc = struct(pot.scgpotc{i});
    p(i) = potc.p;
    if ~isempty(n1i)
        A1(:,i) = potc.A(n1i);
    end 
    if ~isempty(n2i)
        A2(:,i) = potc.A(n2i);
    end
    C11(:,:,i) = potc.C(n1i, n1i);
    C12(:,:,i) = potc.C(n1i, n2i);
    C21(:,:,i) = potc.C(n2i, n1i);
    C22(:,:,i) = potc.C(n2i, n2i);
    if isequal(0, C22(:,:,i)) & isequal(evidence{node}, A2(:, i))
        indic = i;
    end
end

np = zeros(1,dsize);
nA = zeros(cheadsize, dsize);
nC = zeros(cheadsize, cheadsize, dsize);

if indic
    np(:) = 0;
    np(indic) = p(indic);
    nA = A1;
    nC = C11;
else
    for i=1:dsize
        if isequal(0, C22(:,:,i))
            p(i) = 0;
            nA(:, i) = A1(:, i);
            nC(:,:,i) = C11(:,:,i);
        else
            sq = (y2 - A2(:,i))' * inv(C22(:,:,i)) * (y2 - A2(:,i));
            ex = exp(-0.5*sq);
            %np(i) = p(i) * ex / ( (2 * pi)^(-onodesize/2) * sqrt(det(C22(:,:,i))) );
            np(i) = p(i) * ex / ( (2 * pi)^(onodesize/2) * sqrt(det(C22(:,:,i))) );
            nA(:,i) = A1(:,i) + C12(:,:,i) * inv(C22(:,:,i)) * (y2 - A2(:,i));
            tmp1 = C12(:,:,i) * inv(C22(:,:,i)) * C21(:,:,i);
            nC(:,:,i) = C11(:,:,i) - tmp1;
        end
    end
end 

scpot = cell(1, dsize);
W = zeros(cheadsize,ctailsize);
for i=1:dsize
    scpot{i} = scgcpot(cheadsize, ctailsize, np(i), nA(:,i), W, nC(:,:,i));
end
ns(node) = 0;
newscgpot = scgpot(ddom, cheaddom, ctaildom, ns, scpot);
