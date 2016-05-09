BEGIN{

  # read in the filling tiles
  i=0;
  while (getline < fillingTimes)
  {
    i++;
    Ft[i] = $2;
  }
  close(fillingTimes);
  Ntime = i;
  
  # read in the retention tiles
  i=0;
  while (getline < retentionTimes)
  {
    i++;
    Rt[i] = $1;
  }
  close(retentionTimes);

  # boundary conditions
  Rt[0]=Rt[1]
  Rt[Ntime+1]=Rt[Ntime]

  i=0;
  max=-1;
}

/TITLE/{ i++; }

{
  if(NF==5 || NF==4) {
     a[$4] += $2/Ft[i] * (Rt[i+1]-Rt[i-1]);
     #print Ft[i]" "Rt[i+1]" "Rt[i-1]
  }   
  if($4>max) max=$4;
} 

END{
     for(i=1;i<=max;i++) 
	print i"\t"a[i]; 
}
