#!/usr/bin/perl -w
use File::Copy;
use File::Path;
use File::Find;

# temperature-information
@temperature = (120,180,240,280,300,320); # list of temperatures
@molecule = ("propane"); # list of molecules

@epsilon=(51,51.5,52,52.5,53);
@sigma=(3.975,3.98,3.985);

# create the results directory
mkpath("Results");

open(DATw1, ">Results/Results.dat") || die("Could not open file!");

# loop over frameworks, temperatures, molecules
$index_molecule=0;
foreach (@molecule)
{
  $index_sigma=0;
  foreach (@sigma)
  {
    $index_epsilon=0;
    foreach (@epsilon)
    {
      open(DATw2, ">Results/Results.dat-$molecule[$index_molecule]-$sigma[$index_sigma]-$epsilon[$index_epsilon]") || die("Could not open file!");
      $index_temperature=0;
      foreach (@temperature)
      {
        $dir_molecule="$molecule[$index_molecule]/$temperature[$index_temperature]K/$sigma[$index_sigma]/$epsilon[$index_epsilon]";

        #search for all RASPA output-files
        @files=();
        find (\&search, "$dir_molecule");
        sub search {push (@files,$File::Find::name) if(-f and /output*/);}

        foreach $file (@files)
        {
          $density=`gawk '/Average Density:/{getline; getline; getline; getline; getline; getline; getline; getline; ++c;if(c==1) print \$2" "\$5" "}' '$file'`;

          printf DATw1 "$temperature[$index_temperature] $density $sigma[$index_sigma] $epsilon[$index_epsilon]";
          printf DATw2 "$temperature[$index_temperature] $density\n";
        } 
        $index_temperature=$index_temperature+1;
      }
      close(DATw2);
      $index_epsilon=$index_epsilon+1;
    }
    $index_sigma=$index_sigma+1;
  }
  $index_molecule=$index_molecule+1;
}
close(DATw1);

