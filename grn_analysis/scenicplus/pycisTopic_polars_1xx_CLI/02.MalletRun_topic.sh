mallet_path='~/work/Mallet/bin/mallet'
mallet_corpus_path='~/work/MalletRun/tmp/corpus.mallet'
output_prefix='~/work/MalletRun/MalletRun'

## run topic modeling
pycistopic topic_modeling mallet run -i ${mallet_corpus_path} -o ${output_prefix} -t 15 20 25 30 35 40 45 -p 30 -n 150 -a 50 -A True -e 0.1 -E False -s 555 -m 1200 -b ${mallet_path} -v
