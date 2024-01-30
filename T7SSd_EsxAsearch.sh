#!/bin/bash
echo Analysing 41 T7SSd loci using a DNA hmm profile built with hmmbuild to identify EsxA homologues.

assembly=(WP_162178969.1.gb	WP_278773440.1.gb
WP_006290627.1.gb	WP_204666332.1.gb	WP_281667415.1.gb
WP_006290899.1.gb	WP_249113669.1.gb	WP_282742024.1.gb
WP_006720310.1.gb	WP_270472018.1.gb	WP_289613135.1.gb
WP_007147986.1.gb	WP_270574775.1.gb	WP_291288263.1.gb
WP_013252383.1.gb	WP_273061482.1.gb	WP_291289938.1.gb
WP_036575543.1.gb	WP_273383871.1.gb	WP_297035097.1.gb
WP_040590333.1.gb	WP_276572621.1.gb	WP_297183845.1.gb
WP_052118160.1.gb	WP_276867649.1.gb	WP_302153079.1.gb
WP_058271108.1.gb	WP_277142090.1.gb	WP_303942531.1.gb
WP_070638651.1.gb	WP_277150349.1.gb	WP_307388971.1.gb
WP_102371681.1.gb	WP_277156185.1.gb	WP_307738724.1.gb
WP_102372752.1.gb	WP_277157947.1.gb	WP_308615622.1.gb
WP_130810257.1.gb	WP_277166581.1.gb
WP_140396497.1.gb	WP_277175587.1.gb)

strain=(WP_162178969.1.gb	WP_278773440.1.gb
WP_006290627.1.gb	WP_204666332.1.gb	WP_281667415.1.gb
WP_006290899.1.gb	WP_249113669.1.gb	WP_282742024.1.gb
WP_006720310.1.gb	WP_270472018.1.gb	WP_289613135.1.gb
WP_007147986.1.gb	WP_270574775.1.gb	WP_291288263.1.gb
WP_013252383.1.gb	WP_273061482.1.gb	WP_291289938.1.gb
WP_036575543.1.gb	WP_273383871.1.gb	WP_297035097.1.gb
WP_040590333.1.gb	WP_276572621.1.gb	WP_297183845.1.gb
WP_052118160.1.gb	WP_276867649.1.gb	WP_302153079.1.gb
WP_058271108.1.gb	WP_277142090.1.gb	WP_303942531.1.gb
WP_070638651.1.gb	WP_277150349.1.gb	WP_307388971.1.gb
WP_102371681.1.gb	WP_277156185.1.gb	WP_307738724.1.gb
WP_102372752.1.gb	WP_277157947.1.gb	WP_308615622.1.gb
WP_130810257.1.gb	WP_277166581.1.gb
WP_140396497.1.gb	WP_277175587.1.gb)


#locate EsxA gene sequences in the current strain
hmmbuild tsdA_search.hmm esxA_aln.fasta

for index in $(seq 0 40);

do 

#index the current strain for manipulation with HMMER's accessory tool Easel
esl-sfetch --index ${assembly[index]}

#locate tsdA gene sequences in the current strain
nhmmer --tblout ${strain[index]}_tsdA.tblout TsdA_search.hmm ${assembly[index]}

#extract all tsdA-like sequences with an E-value < 5.0e-2 (0.05)
grep -v "^#" ${strain[index]}_tsdA.tblout | awk '{if ($13 < 5.0e-2){print $1"/"$9"-"$10, $9, $10, $1}}' | esl-sfetch -Cf ${assembly[index]} - > ${strain[index]}_tsdA.fa

#add in turn each tsdA-like sequence to a fasta file
cat ${strain[index]}_tsdA.fa >> tsdA_all.fa

done
