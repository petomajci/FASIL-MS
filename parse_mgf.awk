BEGIN {


# read insertions
i=0;
while (getline < spectraIDsFile)
{
  i++;
  spectraID[i] = $1;
}
close(spectraIDsFile);
Nspectra = i;

i=0;
while (getline < Bions)
{
  i++;
  Bion_mz[i] = $1;
}
close(Bions);
Nbions=i;

a=100000000000
i=1;
}

{
#print spectraID[i];
 if(i<=Nspectra && $0 == spectraID[i]) { 
   a=NR
   print spectraID[i]; i++; 
 }

 if ($0=="END IONS")  {
    if (a!=100000000000){
      for(j=1;j<=Nbions;j++){
	if (found[j]==0)  print Bion_mz[j]" 0 ?\t"j"\t"
	found[j]=0;
      }
    }
    a=100000000000
 }

 if (NR==a+4) print
 if (NR>=a+5) {
    for(j=1;j<=Nbions;j++)
	if (Bion_mz[j]-0.015 <= $1 && $1 <= Bion_mz[j]+0.015){
	    found[j]=1;
	    print $0"\t"j"\t"sqrt((Bion_mz[j]-$1)*(Bion_mz[j]-$1)); 
        }
 }

}
