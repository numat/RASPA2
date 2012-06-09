#!/usr/bin/perl -w
use File::Copy;
use File::Path;

# framework-information
@framework = ("MIL-100-primitive", "MIL-101-primitive");  # list of structures
@HeliumVoidFraction = ("0.712419", "0.808882");  # list of helium-voidfractions for the structures
@UnitCells = ("1 1 1","1 1 1"); # list of sizes of the unit cells for the structures
@Forcefield=("GenericMOFs","GenericMOFs");
@RemoveAtomNumberCodeFromLabel=("yes","yes");

# temperature-information
@temperature = (433.0); # list of temperatures
@molecule = ("hexane","2-methylpentane","3-methylpentane","22-dimethylbutane","23-dimethylbutane"); # list of molecules
@idealgas = ([0.00874526,0.0477596,0.0535328,0.226576,0.0873119]); # list of IG Rosenbluth weights for each temperature

# pressure-information
$pressure_start = 1e-5;
$pressure_end = 1e5;
$number_of_pressure_points = 11; # 11 point equally spaced in log-scale
$pressure_unit="Pa";  # "bar", "kPa", or "Pa"
$pressure_scale="log"; # "log" or "linear"

# simulation-information
$SimulationType="MonteCarlo";
$NumberOfCycles="150000";
$NumberOfInitializationCycles="50000";
$PrintEvery="5000";
$RestartFile="no";

# system and queuing information
$divide_into_batches="no"; # combine serial run in larger blocks
$batches = 8; # combine into an 8-core job
$queue = "mof1"; # the queue type
$job_name = "Mixture"; # name of the job
@file_list = (); # list of files copied to all the directories

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

if($pressure_scale eq "log")
{
  for($i = 0; $i < $number_of_pressure_points; $i++)
  {
    $pressure[$i]=10**(((log10($pressure_end)-log10($pressure_start))*($i/($number_of_pressure_points-1)))+log10($pressure_start));
  }
}
else
{
  for($i = 0; $i < $number_of_pressure_points; $i++)
  {
    $pressure[$i]=($pressure_end-$pressure_start)*($i/($number_of_pressure_points-1))+$pressure_start;
  }
}

# get cluster-name
chomp($cluster = `hostname -s`);

# empty 'submit-file'
open(DATw4, ">submit") || die("Could not open file!");
close(DATw4);

$index_batches=0;
$count=0;

$index_framework=0;
foreach (@framework)
{
  $dir_fr=$framework[$index_framework];
  mkpath($dir_fr);

  $index_temperature=0;
  foreach (@temperature)
  {
    # concat the molecule names with "-" between them
    $molecule_list = join("-",@molecule);

    $dir_temp="$dir_fr/$temperature[$index_temperature]K/$molecule_list";
    mkpath($dir_temp);

    $index_pressure=0;
    foreach (@pressure)
    {
      $dir_press="$dir_temp/$pressure[$index_pressure]$pressure_unit";
      mkpath($dir_press);

      open(DATw1, ">$dir_press/run") || die("Could not open file!");
      printf DATw1 "#! /bin/sh -f\n";
      printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
      printf DATw1 "\${RASPA_DIR}/bin/simulate \$1";
      close(DATw1);
      chmod 0755, "$dir_press/run";

      # copy files to the dir
      for($i = 0; $i < scalar @file_list; $i++)
      {
        copy("$file_list[$i]","$dir_press/$file_list[$i]") or die "Copy failed: $!";
      }

      if(($cluster eq "login3") || ($cluster eq "login4"))
      {
        open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
        printf DATw1 "#!/bin/bash\n";
        printf DATw1 "#PBS -S /bin/bash\n";
        printf DATw1 "#PBS -l nodes=1\n";
        printf DATw1 "#PBS -N $job_name\n";
        printf DATw1 "#PBS -l walltime=120:00:00\n";
        printf DATw1 "#PBS -o pbs.out\n";
        printf DATw1 "#PBS -e pbs.err\n";
        printf DATw1 "#PBS -r n\n";
        printf DATw1 "#PBS -V\n";
        printf DATw1 "#PBS -mba\n";
        printf DATw1 "\n";
        printf DATw1 "cd \${PBS_O_WORKDIR}\n";
        printf DATw1 "\n";
        printf DATw1 "echo \$JOB_ID > jobid\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin//simulate \$1\n";
        close(DATw1);
        chmod 0755, "$dir_press/bsub.job";
      }
      elsif($cluster eq "carbon")
      {
        open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
        printf DATw1 "#!/bin/bash\n";
        printf DATw1 "# script for Grid Engine\n";
        printf DATw1 "#\$ -S /bin/bash\n";
        printf DATw1 "#\$ -q $queue\n";
        printf DATw1 "#\$ -N $job_name\n";
        printf DATw1 "#\$ -V\n";
        printf DATw1 "#\$ -cwd\n";
        printf DATw1 "\n";
        printf DATw1 "echo \$JOB_ID > jobid\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
        close(DATw1);
        chmod 0755, "$dir_press/bsub.job";
      }
      elsif($cluster eq "kraken")
      {
        open(DATw1, ">$dir_press/bsub.job") || die("Could not open file!");
        printf DATw1 "#!/bin/bash\n";
        printf DATw1 "# script for Grid Engine\n";
        printf DATw1 "#\$ -S /bin/bash\n";
        printf DATw1 "#\$ -q $queue\n";
        printf DATw1 "#\$ -N $job_name\n";
        printf DATw1 "#\$ -V\n";
        printf DATw1 "#\$ -cwd\n";
        printf DATw1 "\n";
        printf DATw1 "# set path for intel compiler libraries\n";
        printf DATw1 "source /opt/intel/cce/10.1.018/bin/iccvars.sh\n";
        printf DATw1 "\n";
        printf DATw1 "echo \$JOB_ID > jobid\n";
        printf DATw1 "export RASPA_DIR=\${HOME}/RASPA/simulations/\n";
        printf DATw1 "\${RASPA_DIR}/bin/simulate \$1\n";
        close(DATw1);
        chmod 0755, "$dir_press/bsub.job";
      }
      else
      {
        die("Unknown cluster!");
      }

      
      open(DATw3, ">$dir_press/simulation.input") || die("Could not open file!");
      print DATw3 "SimulationType                $SimulationType\n";
      print DATw3 "NumberOfCycles                $NumberOfCycles\n";
      print DATw3 "NumberOfInitializationCycles  $NumberOfInitializationCycles\n";
      print DATw3 "PrintEvery                    $PrintEvery\n";
      print DATw3 "RestartFile                   $RestartFile\n";
      print DATw3 "\n";
      print DATw3 "ContinueAfterCrash no\n";
      print DATw3 "WriteBinaryRestartFileEvery 5000\n";
      print DATw3 "\n";
      print DATw3 "ChargeMethod                  Ewald\n";
      print DATw3 "Forcefield                    $Forcefield[$index_framework]\n";
      print DATw3 "RemoveAtomNumberCodeFromLabel $RemoveAtomNumberCodeFromLabel[$index_framework]\n";
      print DATw3 "\n";
      print DATw3 "Framework           0\n";
      print DATw3 "FrameworkName       $framework[$index_framework]\n";
      print DATw3 "UnitCells           $UnitCells[$index_framework]\n";
      print DATw3 "HeliumVoidFraction  $HeliumVoidFraction[$index_framework]\n";
      print DATw3 "ExternalTemperature $temperature[$index_temperature]\n";
      if($pressure_unit eq "bar") {printf DATw3 "ExternalPressure    %g\n", $pressure[$index_pressure]*1e5}
      elsif($pressure_unit eq "kPa") {printf DATw3 "ExternalPressure    %g\n", $pressure[$index_pressure]*1e3}
      elsif($pressure_unit eq "Pa")  {printf DATw3 "ExternalPressure    %g\n", $pressure[$index_pressure]}
      print DATw3 "\n";

      $index_molecule=0;
      foreach (@molecule)
      {
        print DATw3 "\n";
        print DATw3 "Component $index_molecule MoleculeName              $molecule[$index_molecule]\n";
        print DATw3 "            StartingBead              0\n";
        print DATw3 "            MoleculeDefinition        TraPPE\n";
        print DATw3 "            IdealGasRosenbluthWeight  $idealgas[$index_temperature][$index_molecule]\n";
        print DATw3 "            FugacityCoefficient       1.0\n";
        print DATw3 "            TranslationProbability    1.0\n";
        print DATw3 "            RotationProbability       1.0\n";
        print DATw3 "            ReinsertionProbability    1.0\n";
        print DATw3 "            CBMCProbability           1.0\n";
        print DATw3 "            IdentityChangeProbability 1.0\n";
        print DATw3 "              NumberOfIdentityChanges ",scalar @molecule,"\n";
        print DATw3 "              IdentityChangesList    ";
        for($i = 0; $i < scalar @molecule; $i++) {print DATw3 " $i"};
        print DATw3 "\n";
        print DATw3 "            SwapProbability           1.0\n";
        print DATw3 "            CreateNumberOfMolecules   0\n";
        $index_molecule=$index_molecule+1;
      }
      close(DATw3);

      if($divide_into_batches eq "yes")
      {
        if($count%$batches==0)
        {
          $index_batches=$index_batches+1;

          if(($cluster eq "login3") || ($cluster eq "login4"))
          {
            open(DATw4, ">submit_$index_batches") || die("Could not open file!");
            print DATw4 "\#\!/bin/bash\n";
            print DATw4 "#PBS -S /bin/bash\n";
            print DATw4 "#PBS -l nodes=1\n";
            print DATw4 "#PBS -N $job_name\n";
            print DATw4 "#PBS -l walltime=120:00:00\n";
            print DATw4 "#PBS -o pbs.out\n";
            print DATw4 "#PBS -e pbs.err\n";
            print DATw4 "#PBS -r n\n";
            print DATw4 "#PBS -V\n";
            print DATw4 "#PBS -mba\n";
            print DATw4 "\n";
            print DATw4 "cd \$\{PBS_O_WORKDIR\}\n";
            print DATw4 "\n";
            close(DATw4);
            chmod 0755, "submit_$index_batches"
          }
          elsif($cluster eq "carbon")
          {
            open(DATw4, ">submit_$index_batches") || die("Could not open file!");
            printf DATw4 "#!/bin/bash\n";
            printf DATw4 "# script for Grid Engine\n";
            printf DATw4 "#\$ -S /bin/bash\n";
            printf DATw4 "#\$ -q $queue\n";
            printf DATw4 "#\$ -N $job_name\n";
            printf DATw4 "#\$ -V\n";
            printf DATw4 "#\$ -cwd\n";
            printf DATw4 "#\$ -pe orte $batches\n";
            printf DATw4 "\n";
            close(DATw4);
            chmod 0755, "$dir_press/bsub.job";
          }
          elsif($cluster eq "kraken")
          {
            open(DATw4, ">submit_$index_batches") || die("Could not open file!");
            printf DATw4 "#!/bin/bash\n";
            printf DATw4 "# script for Grid Engine\n";
            printf DATw4 "#\$ -S /bin/bash\n";
            printf DATw4 "#\$ -q $queue\n";
            printf DATw4 "#\$ -N $job_name\n";
            printf DATw4 "#\$ -V\n";
            printf DATw4 "#\$ -cwd\n";
            printf DATw4 "#\$ -pe orte $batches\n";
            printf DATw4 "\n";
            printf DATw4 "# set path for intel compiler libraries\n";
            printf DATw4 "source /opt/intel/cce/10.1.018/bin/iccvars.sh\n";
            printf DATw4 "\n";
            close(DATw4);
            chmod 0755, "$dir_press/bsub.job";
          }
          else
          {
            die("Unknown cluster!");
          }

          open(DATw5, ">>submit") || die("Could not open file!");
          print DATw5 "qsub ./submit_$index_batches\n";
          close(DATw5);
        }

        open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
        print DATw5 "cd \"$dir_press\"\n";
        print DATw5 ".\/run > out 2\>\&1\&\n";
        print DATw5 "cd -\n";

        if(($count+1)%$batches==0) {print DATw5 "wait\n";}
        close(DATw5);

        $count=$count+1;
      }
      else
      {
        open(DATw4, ">>submit") || die("Could not open file!");
        print DATw4 "cd \"$dir_press\"\n";
        print DATw4 "qsub bsub.job\n";
        print DATw4 "cd -\n";
        close(DATw4);
      }


      $index_pressure=$index_pressure+1;
    }
    $index_temperature=$index_temperature+1;
  }
  $index_framework=$index_framework+1;
}

# put a "wait" on the last submit which (probably) is not a full batch
if($divide_into_batches eq "yes")
{
  open(DATw5, ">>submit_$index_batches") || die("Could not open file!");
  if(($count)%$batches!=0) {print DATw5 "wait\n";}
  close(DATw5);
}

chmod 0755, "submit"

