# Run with a single K (e.g. K=11) not a range
#./Genetic\ Analysis/ADMIXTURE/admixture --cv PLINK/gl_plink.bed 1 > ./Genetic\ Analysis/ADMIXTURE/outputs/1.out

# Run with K=1 to K=11 10 runs each
for k in {1..11}; do for r in {1..10};
do
  ./Genetic\ Analysis/ADMIXTURE/admixture -s ${RANDOM} --cv PLINK/gl_plink.bed $k > ./Genetic\ Analysis/ADMIXTURE/outputs/logK${k}r${r}.out
  mv gl_plink.${k}.Q gl_plink.${k}r${r}.Q 
done; done

# Check cross-validation error - lowest is best
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//'

# pong to visualise
pong -m ./outputs/pong_filemap.tdv -i ./outputs/ind2pop.txt

#pong is freely available and can be installed using the Python
#package management system pip. pong’s source code is available at https://github.com/abehr/pong
