#!/bin/bash -l
#########################################################################################################################
########## SkimSeq_Method, average read count per chromosome arm based on provided centromeric position ################# 
#SBATCH --job-name=Skim-Seq-Manuscript
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=5G    # Memory per core, use --mem= for memory per node
#SBATCH --output="%x_%j.o"
#SBATCH --error="%x_%j.e"

for filename in *.txt; ## the input text files should have at least three columns: read count, chromosome, and Position information
do

echo "Output of file $filename ">>Ratio-Output
echo "chromosome    Short Arm Ratio     Long Arm Ratio  ">>Ratio-Output

echo "Processing chr1A from file $filename"
cat $filename|grep chr1A >>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>213000000) next }1' temp>>shortarmtemp   ## fixed wheat chromosome 1A centromere at 213 MB
awk 'NF>1 {if ($3<=213000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/213" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -213" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr1A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr1B from file $filename"
cat $filename|grep chr1B >>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>240000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=240000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/240" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -240" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr1B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr1D from file $filename"
cat $filename|grep chr1D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>170000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=170000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/170" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -170" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr1D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr2A from file $filename"
cat $filename|grep chr2A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>340000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=340000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/340" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -340" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr2A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr2B from file $filename"
cat $filename|grep chr2B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>348000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=348000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/348" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -348" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr2B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr2D from file $filename"
cat $filename|grep chr2D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>268000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=268000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/268" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -268" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr2D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr3A from file $filename"
cat $filename|grep chr3A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>317000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=317000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/317" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -317" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr3A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr3B from file $filename"
cat $filename|grep chr3B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>348000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=348000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/348" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -348" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr3B         $sarmratio                      $larmratio">>Ratio-Output
echo "Processing chr3D from file $filename"
cat $filename|grep chr3D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>240000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=240000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/240" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -240" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr3D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr4A from file $filename"
cat $filename|grep chr4A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>315000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=315000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/315" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -315" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr4A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr4B from file $filename"
cat $filename|grep chr4B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>318000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=318000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/318" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -318" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr4B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr4D from file $filename"
cat $filename|grep chr4D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>185000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=185000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/185" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -185" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr4D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr5A from file $filename"
cat $filename|grep chr5A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>250000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=250000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/250" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -250" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr5A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr5B from file $filename"
cat $filename|grep chr5B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>200000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=200000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/200" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -200" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr5B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr5D from file $filename"
cat $filename|grep chr5D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>186000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=186000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/186" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -186" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr5D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr6A from file $filename"
cat $filename|grep chr6A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>291000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=291000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/291" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -291" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr6A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr6B from file $filename"
cat $filename|grep chr6B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>325000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=325000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/325" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -325" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr6B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr6D from file $filename"
cat $filename|grep chr6D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>214000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=214000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/214" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -214" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr6D         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr7A from file $filename"
cat $filename|grep chr7A>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>362000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=362000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/362" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -362" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr7A         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr7B from file $filename"
cat $filename|grep chr7B>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>294000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=294000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/294" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -294" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr7B         $sarmratio                      $larmratio">>Ratio-Output

echo "Processing chr7D from file $filename"
cat $filename|grep chr7D>>temp
total="$(sum=0;awk '{sum+=$1;}END{print sum;}' temp)"
awk 'NF>1 {if ($3>338000000) next }1' temp>>shortarmtemp
awk 'NF>1 {if ($3<=338000000) next }1' temp>>longarmtemp
rm temp;
sarmsum=0;larmsum=0;
sarmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' shortarmtemp)"
sarmratio=$(echo "scale=2;$sarmsum/338" | bc -l )
rm shortarmtemp;
larmsum="$(sum=0;awk '{sum+=$1;}END{print sum;}' longarmtemp)"
largestvalue="$(tail -1 longarmtemp | awk '{print $3}')"
largestvalue=$(echo "scale=2;$largestvalue/1000000 -338" | bc -l )
larmratio=$(echo "scale=2;$larmsum/$largestvalue" | bc -l )
rm longarmtemp;
echo "chr7D         $sarmratio                      $larmratio">>Ratio-Output

cat Ratio-Output > ${filename}.chrcount.per.Mb.txt
done
