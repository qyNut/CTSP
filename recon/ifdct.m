function x = ifdct(C, fdctPara)
    
    is_real = fdctPara.ifdct_is_real;
    M = fdctPara.M;
    N = fdctPara.N;
    x = ifdct_wrapping(C, is_real, M, N);

end