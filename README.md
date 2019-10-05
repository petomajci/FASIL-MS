
FASIL-MS
A set of analysis scripts to relatively quantify all possible combinations of heavy-light acetylations in a peptide containing multiple lysines for data acquired by FASIL-MS method.

Author: Peter Májek
contact: pmajek@cemm.at

Here we provide a set of scripts to solve for relative abundance of possible heavy-light acetylation combinations of a peptide. The programs can be executed on a linux platform with installed R (www.r-project.org). The package mgcv has to be installed under R. 

Calculation is executed as

./solve.sh <MGF_file> <begin_retention_time> <end_retention_time>

The programs will integrate intensities of all relevant b- and y- ions in the retention time window (specified on the command line in minutes).

The relevant peptide sequence shall be specified in the file peptide.sequence and masses of all y- and b- ions shall be specified in the file peptide.modifications. Please inspect the later file for the details of the format. The masses of y- and b-ions in the peptide.modifications file have been obtained from this online fragment ion calculator: http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html. If you desire not to use some of the b- or y- ions for the analysis (for example due to interfering internal fragments) then you still have to keep those ions in the peptide.modifications, just please change the mass of the ions you desire to exclude from the analysis to some dummy mass not present in your spectra. In the provided file this was done for b-3 and b-4 ions by changing their fragment masses like 285.15676 => 11285.15676. 

The last piece of input is the filling time of each MS2 spectra this shall be provided in the file spectra_filling_times.txt (examine the example file for the details of the format). We have used Proteome discoverer to extract the filling times from raw files.

The generated result file will list relative abundances of all heavy light acetylation combinations specified in the file peptide.modifications, in the same order as given in that file. We used Proteowizzard to generate mgf files from thermo raw files.

files: M495-E05-C643-P7617-1.mgf spectra_filling_times.txt M495-E05-C643-P7617-1-result.txt are provided as examples for calculation of H4 peptide acetylation combinations (specified in peptide.sequence, peptide.modifications)

The command ./solve.sh M495-E05-C643-P7617-1.mgf 38.50 39.80 shall re-generate the result file M495-E05-C643-P7617-1-result.txt.

You would need to change content of peptide.sequence, peptide.modifications, and spectra_filling_times.txt if you want to analyze your peptide of interest.

We also provide an auxiliary script check.pl that can be used to check for theoretically possible internal fragmet ions having the same m/z as one of the y or b ions. It is executed as "perl check.pl" and it check the peptide sequence listed in the file peptide.sequnce.
