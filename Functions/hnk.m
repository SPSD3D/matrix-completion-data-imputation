function  x = nhk(b,l)
%This Function is used to Create Hankel Type Data Matrix
%x is a given data
%l represent window size
    n = length(b);
    m = n-l+1;
    x = zeros(m,l);
    for i=1:m
        x(i,:)=b(i:i+l-1);
    end;
