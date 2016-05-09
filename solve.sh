#!/bin/bash

INPUT=$1

startTime=$(bc <<< "scale = 3; $2 * 60 ")
endTime=$(bc <<< "scale = 3; $3 * 60 ")

PDpsms=spectra_filling_times.txt

sequence=$(cat peptide.sequence)
seqL=${#sequence}
 L=$((${#INPUT}-4))
MGFname=${INPUT:0:L}

CODE=$MGFname

MODIFfile=peptide.modifications 

N=$(echo $sequence | grep -o "K" | wc -l)
Nstates=$(($N*$N))

MGF=$MGFname.mgf

grep -P "RTINSECONDS=|TITLE=" $MGF | awk '{if(NR%2==0) print; else printf("%s\t",$0)}' | sed -e "s/RTINSECONDS=//g" | awk -v s=$startTime -v e=$endTime '{ if(s<=$2 && $2<=e) print $1}'  > $CODE.psmIDs
cat $CODE.psmIDs | sed -e "s/TITLE=//g" -e "s/dta.*/dta/g" > $CODE.tmp
grep -P "RTINSECONDS=|TITLE=" $MGF | awk '{if(NR%2==0) print; else printf("%s\t",$0)}' | sed -e "s/RTINSECONDS=//g" | awk -v s=$startTime -v e=$endTime '{ if(s<=$2 && $2<=e) print $2}' > $CODE.psmRTs


grep $MGFname.raw $PDpsms | sed -e "s/\"//g" -e "s/.raw//g"| awk '{printf("%s.%s.%s.%s.dta\t%s\n",$4,$3,$3,$2,$1)}' | sort | uniq > all.fillingTimes
# we have to check if we have filling time in all.fillingTimes file as for some high charge states (8+) it is missing, I just use filling time of 200 then
for f in $(cat $CODE.tmp); do A=$(grep $f all.fillingTimes); if [ "$A" == "" ]; then A="$f\t200"; fi; echo -e "$A"; done > $CODE.psmTimes

for((i=1;i<=$seqL;i++)); do awk -v i=$i '{if($2==i) print }' $MODIFfile | awk '{a[NR]=$3} END{N=NR; for(i=1;i<=N;i++) { for(j=1;j<=N;j++) if(a[i]==a[j]) printf("1 "); else printf("0 "); print ""; }    }'; done  | sed -e "s/ $//g" > B_ions.matrix

for((i=1;i<=$seqL;i++)); do awk -v i=$i '{if($2==i) print }' $MODIFfile | awk '{a[NR]=$4} END{N=NR; for(i=1;i<=N;i++) { for(j=1;j<=N;j++) if(a[i]==a[j]) printf("1 "); else printf("0 "); print ""; }    }'; done | sed -e "s/ $//g" > Y_ions.matrix


for((i=1;i<=$seqL;i++)); do awk -v i=$i '{if($2==i) print }' $MODIFfile | awk '{ print $3 }'; done | sed -e 's/ $//g' > B_ions.mz
for((i=1;i<=$seqL;i++)); do awk -v i=$i '{if($2==i) print }' $MODIFfile | awk '{ print $4 }'; done | sed -e 's/ $//g'> Y_ions.mz


# B-ions problem

sort B_ions.mz | uniq > uniq_B_ions.mz
awk -v spectraIDsFile=$CODE.psmIDs -v Bions=uniq_B_ions.mz -f parse_mgf.awk $MGF > aaa
awk -v fillingTimes=$CODE.psmTimes -v retentionTimes=$CODE.psmRTs -f merge_intensities.awk aaa > bbb
paste uniq_B_ions.mz bbb > ccc
for f in $(cat B_ions.mz); do grep -P "^$f" ccc; done | sed -e 's/\t$/\t0.1/g' > B_ions.mz.intensities
for((i=1;i<=$seqL;i++)); do A=$(($i*$Nstates)); head -n $A B_ions.mz.intensities | tail -n $Nstates | sort | uniq | awk '{a+=$3} END{print a}'; done > weightsB.txt

# Y-ions problem
sort Y_ions.mz | uniq > uniq_Y_ions.mz
awk -v spectraIDsFile=$CODE.psmIDs -v Bions=uniq_Y_ions.mz -f parse_mgf.awk $MGF > aaa
awk -v fillingTimes=$CODE.psmTimes -v retentionTimes=$CODE.psmRTs -f merge_intensities.awk aaa > bbb
paste uniq_Y_ions.mz bbb > ccc
for f in $(cat Y_ions.mz); do grep -P "^$f" ccc; done | sed -e 's/\t$/\t0.1/g' > Y_ions.mz.intensities
for((i=1;i<=$seqL;i++)); do A=$(($i*$Nstates)); head -n $A Y_ions.mz.intensities | tail -n $Nstates | sort | uniq | awk '{a+=$3} END{print a}'; done > weightsY.txt

Rscript solveOptimization.R
mv result.txt $CODE-result.txt

# clean

rm $CODE.psmIDs $CODE.tmp $CODE.psmRTs $CODE.psmTimes all.fillingTimes B_ions.matrix Y_ions.matrix B_ions.mz Y_ions.mz uniq_B_ions.mz B_ions.mz.intensities weightsB.txt uniq_Y_ions.mz aaa bbb ccc Y_ions.mz.intensities weightsY.txt 
