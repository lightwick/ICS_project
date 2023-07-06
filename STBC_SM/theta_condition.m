% In 16-QAM, n<=6 needs to hold true for 'theta_k=(k-1)*pi/(2*n)'

for i=3:20
    c = find_c(i);
    a = floor(i/2);
    n = ceil(c/a);

    disp('for n_t:'+string(i)+' n='+string(n));
end

function k = find_c(n_t)
    k = n_t*(n_t-1)/2;
    while bitand(k,k-1)
        k=k-1;
    end
end