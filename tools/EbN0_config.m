assert(exist(EbN0_dB, 'var'), 'EbN0_dB must exist in the workspace')

EbN0 = db2pow(EbN0_dB);

EsN0 = EbN0 * log2(M);
EsN0_db = pow2db(EsN0);