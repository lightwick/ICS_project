function EbN0_dB = EsN0_dB2EbN0_dB(M, EsN0_dB)
    EsN0 = db2pow(EsN0_dB);
    EbN0 = EsN0 / log2(M);
    EbN0_dB = pow2db(EbN0);
end