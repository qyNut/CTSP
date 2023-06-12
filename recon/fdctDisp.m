function fdctDisp(x, fdctPara)

    C = fdct(x, fdctPara);
    coeffdisp_original = fdct_wrapping_dispcoefMN(C,fdctPara.M,fdctPara.N);
    figure;
    imshow(abs(coeffdisp_original));
    
end

  