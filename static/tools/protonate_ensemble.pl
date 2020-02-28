#!/usr/bin/perl

print STDERR"
##########################################################
# protonate_ensemble.pl [-c<command>] < in.pdb > out.pdb #
#                                                        #
# adds hydrogen atoms to protein structures in an        #
# ensemble. Steps performed:                             #
# 1) separation of the ensemble to single PDB structures #
#    based on MODEL/ENDMDL records                       #
# 2) invoking the pdb2gmx program of GROMACS to add      #
#    hydrogens (oplsaa force field, -ignh flag, you      #
#    can modify this with option -c or in the code)      #
# 3) concatenating the protonated structures to generate #
#    the output ensemble                                 #
#                                                        #
# Contact: gaspari.zoltan@itk.ppke.hu                    #
##########################################################

";

use Getopt::Std;

# You can rewrite the default GROMACS command here. "model_in.pdb" represents
# the input pdb file from the original ensemble, "model_out.pdb" 
# stands for the output single-model pdb file. Do _not_ change their names here.

#$opt_c="pdb2gmx_mumo -f model_in.pdb -o model_out.pdb -ff oplsaa -ignh -water none";

$opt_c="pdb2gmx_mumo -f model_in.pdb -o model_out.gro -ff oplsaa -ignh -water none >&/dev/null  \n editconf_mumo -f model_out.gro -o model_out.pdb";



getopts('chv') || die "Unknown/wrong options, please try again\n";

if ($opt_h){die "Use -v for verbose output and -h to display this message.\n";}

print STDERR "Will use the following GROMACS command:\n$opt_c\n" if $opt_v;

$firstremark=1;
$modelnum=0;
while(<>){
    chomp;
    $inline=$_;
    if ($inline=~/^MODEL/){
	&checkfirstremark;
	$modelnum++;
	print "$_\n";
	$atomlines="";
    }
    elsif ($inline=~/^ATOM/){
	&checkfirstremark;
	$atomlines.="$_\n";
    }
    elsif ($inline=~/^ENDMDL/){
	open (IP,">model_in.pdb") || die "Cannot create file model_in.pdb\n";
	print IP $atomlines;
	close(IP);
        # deleting backup files, if any, in order to not overwhelm gromacs
        # backuping (may boil out after 99 backups, so better avoid this)
	system("rm \\#\*") if (-e "\#model_out.pdb.1\#");
        #executing the protonation command, 
 	system("$opt_c >& /dev/null");
	open (OP,"model_out.pdb") || die "Cannot open file model_out.pdb\n";
	while(<OP>){
	    chomp;
	    print "$_\n" if ($_=~/^ATOM/);
	}
	close(OP);
	print "ENDMDL\n";
    }
    elsif ($_=~/^REMARK/){
	&checkfirstremark;
	print "$_\n";
    }
    else{
	print "$inline\n" unless ($inline=~/^ANISOU|^HETATM/);
    }
}

print STDERR "

Found $modelnum models in the input file.
Please check your output manually, especially the chain identifier(s), 
if any were present in the input.
You might also want to convert the nomenclature of the
hydrogen atoms with pdbfile2bmrbnomenclature.pl 
as a next step.

";


sub checkfirstremark{

# Put this line as REMARK before the first REMARK, MODEL or ATOM
# record, whichever present in the input.
    if ($firstremark){
      print "REMARK       PROTONATED WITH THE SCRIPT protonate_ensemble.pl\n";
      $firstremark=0;
    }
}
