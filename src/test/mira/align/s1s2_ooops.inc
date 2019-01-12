
//bla.acquireSequences((const char*) seq1, strlen(seq1),
//             (const char*) seq2, strlen(seq2),
//             1, 2, 1, -1,
//             true, 0);

char seq1[]="AGCAGCGACTGCGCGCACGTATGACTGCATGCTGCGTCTGATGCGCTCTCGCATGTGACGACAGCGCTGCGCTGCGCTCTCTGCGTGAGATACGAGTGTGTCACTGTGCGACTGACGTGCGCGCTGCATGTGCTCTCGACACTGCGTACGTACTGTATCGCANTGACGCGATATCGAGCTGTATCGTGACGCGACTGTACATAGACGAGATATGNACATGCACGCATAGATGTCGCGAGTCGCTGAGATGAGTATCGTGCATCTACTGCNNNNNNNNNNNNNNGCGATCGTCTGACGATGCTCGTGCAGCGATCGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGATGCTGATCTGTATCGCGTCTATCTGTCATCACGCGAGACATGTG";

char seq2[]="GTGTGCATCTGCATCAGCGTATATCACAGCGACATACGTAGTACGATCTGCTCAGTCTGCGTGTCGCTCAGCAGAGTACACTGCTGATGCATCAGCTGTGCTGAGTGATGACTGTCTGATGAGTGCGTCATCTCTGTCAGCGTGCTGACGCTGTGTCATCTGTCGCACTCACTGCTGCACACAGCGCTAGCTGACAGCTCAGCAGAGCAGTCGACGTAGTATCATAGTCGTCTCTCGTCGCATAGCTATGTAGAGTAGCAGCGAGCATCGATCTATCTCAGCGTGCGCTGTGTATCGTCAGTCAGTATAGATCGTACACTGTCAGACGCGTATGCATGCATCTGTCGCGTGACACTAGTATGACTCTCGTCTCGTGCGAGTCATCGTACGCGCGTCGCAGTGTCTGTGTGTATGACGTCACACTGAGTGCGAGCTGTATGCGCGTGTGTAGCACGACGTGCGTGACGTATGAGAGATGCGTAGCGATGCTGCGCGCTGACGTCGTGAGTCGATGACGTCGATGCATCATCATATATCGCAGCTGCGCGCACGTGAGCGCGAGTGCGAGCAGATGATGCATCTGTGATGAGCATGAGCTCGATAGTCGCTATATCGTGCTCTGTCTACGCAGTCGAGATATCAGATGCGAGAGATCGCGCGCTGAGACTAGCGTGCAGACGAGTATACATGCGCGCGCTAGCTGCGCACGAGACGTACACGCATCACGATGTCGTGCGTCATGTCAGATGCGATGTCGTGAGTACACACGTGCGTGACACATGCGAGTACACGCGTGCGTATGCGTGAGCAGTGCGACGCTATGTCACACTCGACTGTAGCATCATGTCGTACGCACGTGCTGTGACGCACTCATCGCGCTGCGTAGCGCGCGTGTGCTCGATGCAGATCTACGCAGTCGCTGCGTGACTACGTCTGCACAGCTGCATGACGTCAGCACGCTGATGATGAGTGAGCGAGCATGATAGCTCTCAGCGATGCGCACGCTCTATATCAGTACTGCGCATCGCTGATGTCGCACGCGATGCGAGCGCGAGCAGTGATCACTCGCAGCGCGCGTGATGTACACTACTGCGCGCGCAGTCATGCATACAGCAGCGCACGTGCTGCGATGATCTGCTACGTGCGTGTCTGACTGATGACACTCAGCGCTCTATCGCACATCTCGTCAGCGATACGTATGTGATGATAGCATCGCGTGACGCGTACATCGCATGTGCAGCAGCGCGACTGTCGCGCGTACGCAGAGATCTCGTGAGTGCGACGTGCGTACTGCTGATCAGAGCATAGTCATCTGATCTGTGTACTGATGATGCACTGTACAGCATACAGCATCTGTGAGCGCACGCTCGCTGCGCATGCTGCAGCGTAGTCGATGAGCTGATGTCTATATCGCACTACTATAGCAGCGACGCGCGATATGCACGCTGCAGCTGAGATAGCGACGCATATCAGAGATACGTCATGTAGCGATGTGCTGCATAGTGTGCATGATCTCGCTCGAGAGCAGCGCACGATATATGCTATCGCTGCGCATACTCGACGAGAGATGACATATGAGCTACGTGCTACGACGCTCGCGTCGCATCACATCGCGCGTCACACTATCATGACAGCTGCTGAGTGCATGACGCACTGACGTATCGTCAGTAGATGACATACGTGTGTGCTGTGCTATGAGCTGCGCTGCGCGTGCACTGCTACGTATACGACGCGTGATGCGCTGTACGTATCGTGTGCACTGCACTGATATGTCGCTGTGTGCAGCACGCTGCGCATGTGCAGACAGCGATCGTGCTGTCTGACACTGAGCATGACAGCGATCGACGTGCTGCACAGCTGCACAGTGCAGAGCTGCACTGACGCGCTGATGATATGTATGAGCATCAGCTGATAGTAGATCGTGACTGCTGCTCGTCGCGCGCTCGTGAGTGTGCGTCAGCGCTCTACTCTGCAGTGTCAGACGACATCGCTGATGTGATACAGTCTGCTGACAGATGTAGACGTGACGTCTGTACTCGTGAGCTGCTGCGTGCGACTATACGCATACTCGACGACTGATGAGCATACTGTCGCTGCTGAGACTGACAGTAGAGCAGTACTGCGCATGCGCTGTACGACTGTCTACGTAGCGATGTGCATACAGTGCATACGAGCGATCGACTGCGATCGTCGCGCAGAGACAGCATAGTCGTCACGCGTACTCGATAGCAGCACTGTGATCGCTACAGTGATATCAGCGTAGATCACGAGACGTCATGATCGCGCTCTACTGCTGAGCATACGTATGCATGTGCGAGTCTCTGATGCGTATGACGTGTAGAGTCTCGTCTGATGTCGACTGCATCGCGACGATGCGCACTCTCATATCGAGACAGACTGCGATCAGCACTGATACGCTGTGCGTACATCTCTGATATCGATCTGCTGATCTGCTATAGCACTGCGTGACTAGCGATCTGCAGCGATGCGATCGCTGCTGACGCTGCATACTACGACGTAGACGAGCGTGCTGAGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACTCATGCGCGATACACGCGTGCATATCATGACGTGTATCACGCTGTCGATGACGAGCTACGCGCTCTGCGACGTGACTGATCTATGTGACGTACACGCGATAGCGCGCGCTGCGATGCACGTGCAGTGCGATGCTAGATGCGACTACTGATCACGTCATCTCGTCATCGCGCGTATACAGCAGCACTGCGCGATCTGCATCAGCATCGATCTCTCTGACGTGATCGATCTCGCATATCGTGAGCAGCGACTGCGCGCACGTATGACTGCATGCTGCGTCTGATGCGCTCTCGCATGTGACGACAGCGCTGCGCTGCGCTCTCTGCGTGAGATACGAGTGTGTCACTGTGCGACTGACGTGCGCGCTGATGTGCTCTCGACACTGCGTACGTACTGTATCGCATGACGCGATATCGAGCTNNNNNNNNNNNNNNNNNNNGACTGAGATATGACATGCACGCATAGATGTCGCGAGTCGCTGAGATGAGTATCGTGCATNNNNNNNNNNNNNNNTCGAGTGATAGCGATCGTCTGACGATGCTCGCTGCAGCGATCGTCTGTATCTGATGCTGCTGCAGCTCTACTCGCTGTGTGATGCTGATCTGTATCGCGTCTATCTGTCATCACGCGAGACATGTGTGCACGATCATGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCGTATCTGTATGTGAGCGTGCGATCGCTGTGATCGCATGTGCATGCGCATGATCANNNNGACGATATGTCGCGTGATCGCAGCATCTGTATATCGTCATCACTGCGCTGTATCTGTCGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATCACGCGTCATAGCTCATCAGCTGTCTCGCGCGAGTCGTATCGTGCGCATGCGTACTCACACAGCAGCGCACTGTGCTCGTGTCGCAGTCTCACTGACGCTGTACGCTCTGCGATCTGCGCGTCNNNNNNNNNNNNNNNNNNNNNNNNGACACAGTCGTGCGTAGTCTCACGACGCGTCGACGCTGCAGCGCTCATGATACGCACGCTCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
