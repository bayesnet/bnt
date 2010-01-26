function [A1, A2, B1, B2, C11, C12, C21, C22] = partition_matrix_vec_3(A, B, C, n1, n2, bs)

dom = myunion(n1, n2);
n1i = block(find_equiv_posns(n1, dom), bs(dom));
n2i = block(find_equiv_posns(n2, dom), bs(dom));


    A1 = A(n1i);
    A2 = A(n2i);
    if isempty(B)
        B1 = zeros(size(n1i, 2),size(B, 2));
        B2 = zeros(size(n2i, 2),size(B, 2));
    else
        B1 = B(n1i, :);
        B2 = B(n2i, :);
    end
    
    
    C11 = C(n1i, n1i);
    C12 = C(n1i, n2i);
    C21 = C(n2i, n1i);
    C22 = C(n2i, n2i);
