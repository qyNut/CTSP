function C = fdct(x, fdctPara)
    
    is_real = fdctPara.fdct_is_real;
    finest = fdctPara.fdct_finest;
    nbscales = fdctPara.fdct_nbscales;
    nbangles_coarse = fdctPara.fdct_nbangles_coarse;

    C = fdct_wrapping(x, is_real, finest, nbscales, nbangles_coarse);

end