assert(exist('EsN0_dB', 'var'), 'EsN0_dB must exist in the workspace')
EsN0 = db2pow(EsN0_dB);

EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);