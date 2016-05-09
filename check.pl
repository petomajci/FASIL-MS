#!/usr/bin/perl
use strict;
use warnings;

my %AAmass;

   $AAmass{'A'} = 89.05;
   $AAmass{'R'} = 174.11;
   $AAmass{'N'} = 132.05;
   $AAmass{'D'} = 133.03;
   $AAmass{'C'} = 121.02;

   $AAmass{'Q'} = 146.07;
   $AAmass{'E'} = 147.05;
   $AAmass{'G'} = 75.07;
   $AAmass{'H'} = 155.07;
   $AAmass{'I'} = 131.09;

   $AAmass{'L'} = 131.09;
   $AAmass{'K'} = 146.11 + 42.01;
   $AAmass{'M'} = 149.06;
   $AAmass{'F'} = 165.08;
   $AAmass{'P'} = 115.06;

   $AAmass{'S'} = 105.04;
   $AAmass{'T'} = 119.06;
   $AAmass{'W'} = 204.09;
   $AAmass{'Y'} = 181.07;
   $AAmass{'V'} = 117.08;


sub sameMass{
  my $tol = 0.05;
  my $tol2 = $tol * $tol;
  my ($M1, $M2) = @_;
  if (($M1 - $M2)*($M1 - $M2) <= $tol2){
    return 1;
  }

 return 0;
}



my $filename = 'peptide.sequence';
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
while (my $peptide = <$fh>) {
  chomp $peptide;

  my $length = length $peptide;

  # prepare b-ion masses array  
  my @BionMass = ();

  # and count of Lysines array
  my @BionLysines = ();

  $BionMass[0] = 0;
  $BionLysines[0] = 0;
  for(my $i=1;$i<=$length-1;$i++){
    $BionMass[$i] = $BionMass[$i-1] + $AAmass{substr($peptide, $i-1, 1)};
    $BionLysines[$i] = $BionLysines[$i-1];
    if (substr($peptide, $i-1, 1) eq 'K'){ $BionLysines[$i] = $BionLysines[$i-1] + 1;
    }
    #print "b-ion mass $BionMass[$i] ",substr($peptide, $i-1, 1),"\n";
  }

  my @conflicts = ();
  for(my $c=0;$c<=$length-1;$c++){
      $conflicts[$c] = 0;
  }

  for(my $i=2;$i<=$length;$i++){
    my $internalFragmentMass = 0; #$AAmass{substr($peptide, $i-1, 1)};
    for(my $j=$i;$j<=$length;$j++){
       $internalFragmentMass += $AAmass{substr($peptide, $j-1, 1)};
       #print substr($peptide, $j-1, 1), "\n";
       
       # compare to all B-ions
         for(my $c=1;$c<=$j-1;$c++){
	   if($BionLysines[$c]>0){

	     for(my $delta=-$BionLysines[$c];$delta<=$BionLysines[$c];$delta++)
	     {
	       if (sameMass($internalFragmentMass, $BionMass[$c]+3.02*$delta)){
	  	  #print "Checking c=$c i=$i j=$j $internalFragmentMass\n";
		  my $foundLysine=0;
	          for(my $ll=1;$ll<$i;$ll++){
		     if (substr($peptide, $ll-1, 1) eq 'K'){ $foundLysine=1; }
  		  }
		  for(my $ll=$c+1;$ll<=$j;$ll++){
		     if (substr($peptide, $ll-1, 1) eq 'K'){ $foundLysine=1; }
	  	  }
	          if ($foundLysine==1 || $delta != 0){
		    $conflicts[$c]++;		
	            #print "Conflict B-ion c=$c i=$i j=$j $internalFragmentMass $BionMass[$c] $delta\n";
		  }
	       }

	     }
           }
         }
     }
  }

  for(my $c=1;$c<=$length-1;$c++){
    if ($conflicts[$c]==1){
        print "There is pottentially $conflicts[$c]", " conflicting internal fragment ions with b-$c ion.\n";
    }
    if ($conflicts[$c]>1){
	print "There are pottentially $conflicts[$c]", " conflicting internal fragment ions with b-$c ion.\n";
    }
  }




  # reverse sequence to check fo potentially conflicting y-ions
  my $revPeptide = reverse $peptide;
  $peptide = $revPeptide;

  # prepare y-ion masses array  
  @BionMass = ();

  # and count of Lysines array
  @BionLysines = ();

  $BionMass[0] = 0;
  $BionLysines[0] = 0;
  for(my $i=1;$i<=$length-1;$i++){
    $BionMass[$i] = $BionMass[$i-1] + $AAmass{substr($peptide, $i-1, 1)};
    $BionLysines[$i] = $BionLysines[$i-1];
    if (substr($peptide, $i-1, 1) eq 'K'){ $BionLysines[$i] = $BionLysines[$i-1] + 1;
    }
    #print "b-ion mass $BionMass[$i] ",substr($peptide, $i-1, 1),"\n";
  }

  @conflicts = ();
  for(my $c=0;$c<=$length-1;$c++){
      $conflicts[$c] = 0;
  }

  for(my $i=2;$i<=$length;$i++){
    my $internalFragmentMass = 0; #$AAmass{substr($peptide, $i-1, 1)};
    for(my $j=$i;$j<=$length;$j++){
       $internalFragmentMass += $AAmass{substr($peptide, $j-1, 1)};
       #print substr($peptide, $j-1, 1), "\n";
       
       # compare to all B-ions
         for(my $c=1;$c<=$j-1;$c++){
	   if($BionLysines[$c]>0){

	     for(my $delta=-$BionLysines[$c];$delta<=$BionLysines[$c];$delta++)
	     {
	       if (sameMass($internalFragmentMass,$BionMass[$c]+3.02*$delta)){
	  	  #print "Checking c=$c i=$i j=$j $internalFragmentMass\n";
		  my $foundLysine=0;
	          for(my $ll=1;$ll<$i;$ll++){
		     if (substr($peptide, $ll-1, 1) eq 'K'){ $foundLysine=1; }
  		  }
		  for(my $ll=$c+1;$ll<=$j;$ll++){
		     if (substr($peptide, $ll-1, 1) eq 'K'){ $foundLysine=1; }
	  	  }
	          if ($foundLysine==1 || $delta != 0){
		    $conflicts[$c]++;		
		    #print "Conflict Y-ion c=$c i=$i j=$j $internalFragmentMass $BionMass[$c] $delta\n";
		  }
	       }

	     }
           }
         }
     }
  }

  for(my $c=1;$c<=$length-1;$c++){
    if ($conflicts[$c]==1){
        print "There is pottentially $conflicts[$c]", " conflicting internal fragment ions with y-$c ion.\n";
    }
    if ($conflicts[$c]>1){
	print "There are pottentially $conflicts[$c]", " conflicting internal fragment ions with y-$c ion.\n";
    }
  }

  # analyze only first line of peptide.sequence file
  exit;
  }
