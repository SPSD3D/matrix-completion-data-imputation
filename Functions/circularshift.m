function C = circularshift( in )

A=in(:);

B=[A(2:end)  ;A(1)];
C=reshape(B,size(in));

end

