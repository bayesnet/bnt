function int = compute_fwd_interface(intra, inter)
% COMPUTE_FWD_INTERFACE Compute nodes with children in the next slice
% function int = compute_fwd_interface(intra, inter)

int = [];
ss = length(intra);
for u=1:ss
  if any(inter(u,:))
    int = [int u];
  end
end
