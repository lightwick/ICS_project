%% Creates all combinations of Nt consisting of number 0~M-1
function Candidates = get_candidates(M, Nt)
    Candidates = zeros(M^Nt, Nt);
     for ii = 0:M^Nt-1
         for jj = 1:Nt
             Candidates(ii+1,jj) = mod(floor(ii/M^(Nt-jj)),M);
         end
     end
end