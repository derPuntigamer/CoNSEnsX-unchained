#!/usr/bin/perl 

###############################################
# invoke the script as:
# perl pdbfile2bmrbnomenclature.pl -h 
# for usage and citing information.


#################################################
# This file contains the minimal versions of 
# two perl classes: Vector3D and PDBStruct
# as well as the code for the actual 
# pdbfile2brmbnomenclature.pl script.
#
# Note that the classes are not intended to be 
# comprehensive, only the routines needed for the
# converter script are included.
#
# The code is provided on an "as is" basis.

##################################################
package Vector3D;

###################################################
## Class for handling 3D points and vectors
## 
## Minimal version to support the script pdbfile2bmrbnomenclature.pl
##
## Version 0.2, 13.02.2002. 
###################################################

sub new {

    my $Vector3D={};
    bless ($Vector3D);
    my $type=shift;
    my $coordx=shift;my $coordy=shift; my $coordz=shift;

    $Vector3D->{COORDX}=$coordx;
    $Vector3D->{COORDY}=$coordy;
    $Vector3D->{COORDZ}=$coordz;

    return ($Vector3D);

}#_sub new

#----------------------------------------------------

sub PrintInfo{

    my $this=shift;
    print "COORDINATES : $this->{COORDX}\t$this->{COORDY}\t$this->{COORDZ}\n";

}#_sub PrintInfo

#----------------------------------------------------

sub Length{

    my $this=shift;
    my $length=sqrt(($this->{COORDX}**2)+($this->{COORDY}**2)+($this->{COORDZ}**2));
    return($length);

}#_sub Length

#----------------------------------------------------

sub Length2{

    my $type=shift;
    my $first=shift;
    return ($first->Length());

}#_sub Length2

#-----------------------------------------------------

# Constructing unity vector from $this

sub Norm{

    my $this=shift;
    my $length=$this->Length();
    my $m_x;  my $m_y;  my $m_z;

    if ($length != 0){
	$m_x=$this->{COORDX}/$length;
	$m_y=$this->{COORDY}/$length;
	$m_z=$this->{COORDZ}/$length;
    }
    else {
	$m_x=0;
	$m_y=0;
	$m_z=0;
    }

    return (new Vector3D($m_x,$m_y,$m_z));

}#_sub Norm

#-----------------------------------------------------

sub Norm2{

    my $type1=shift;
    my $first=shift;
    return ($first->Norm());

}#_sub Norm2

#-------------------------------------------
# Sum vector of $this and another Vector3D object

sub Sum{

    my $this=shift;
    my $second=shift;
    my $m_x=$this->{COORDX}+$second->{COORDX};
    my $m_y=$this->{COORDY}+$second->{COORDY};
    my $m_z=$this->{COORDZ}+$second->{COORDZ};
    return (new Vector3D($m_x,$m_y,$m_z));

}#_sub Sum


#-------------------------------------------

# Sum of two Vector3D objects

sub Sum2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;
    return ($first->Sum($second));

}#_sub Sum2

#-----------------------------------------------------

# Scalar product of $this and another Vector3D object

sub ScalarProduct{

    my $this=shift;
    my $second=shift;
    my $scalarproduct=$this->{COORDX}*$second->{COORDX}+$this->{COORDY}*$second->{COORDY}+$this->{COORDZ}*$second->{COORDZ};

    return ($scalarproduct);

}#_sub ScalarProduct

#-------------------------------------------
# Scalar product of two Vector3D objects

sub ScalarProduct2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;

    return ($first->ScalarProduct($second));

}#_sub ScalarProduct2

#--------------------------------------------

# Vectorial product of $this and another Vector3D object

sub VectorialProduct{

    my $this=shift;
    my $second=shift;
    my $product_x=$this->{COORDY}*$second->{COORDZ}-$this->{COORDZ}*$second->{COORDY};
    my $product_y=$this->{COORDZ}*$second->{COORDX}-$this->{COORDX}*$second->{COORDZ};
    my $product_z=$this->{COORDX}*$second->{COORDY}-$this->{COORDY}*$second->{COORDX};
    return (new Vector3D($product_x,$product_y,$product_z));

}#_sub VectorialProduct


#-------------------------------------------

# Vectorial product of two Vector3D objects

sub VectorialProduct2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;
    return ($first->VectorialProduct($second));

}#_sub VectorialProduct2

#--------------------------------------------


# Difference vector of $this and another Vector3D object

sub Diff{

    my $this=shift;
    my $second=shift;
    my $m_x=$this->{COORDX}-$second->{COORDX};
    my $m_y=$this->{COORDY}-$second->{COORDY};
    my $m_z=$this->{COORDZ}-$second->{COORDZ};
    return (new Vector3D($m_x,$m_y,$m_z));

}#_sub Diff


#-------------------------------------------

# Difference of two Vector3D objects

sub Diff2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;
    return ($first->Diff($second));

}#_sub Diff2

#--------------------------------------------


package PDBStruct;

# Depends on the Vector3D package for operations on atomic
# coordinates
#use Vector3D;

###############################################################
## Class for handling PDB entries
##
## Minimal version to support the script pdbfile2bmrbnomenclature.pl
##
###############################################################

sub new {

    my $PDBObject={};
    bless ($PDBObject);
    my $EntryData=shift;$EntryData=shift;

    $PDBObject->ParseData($EntryData);
    return ($PDBObject);

}#_sub new
#----------------------------------------------------

############################################################
## Function ParseData assigns information provided
## by the PDB file to PDBstruct object variables
############################################################

sub ParseData {


    my $this=shift;
    my $FileData= shift;

## Splitting lines into an array 

    my @FileLines=split(/\n/,$FileData);    
    my $LineRead=0;


    $this->{FILELINES}=\@FileLines;

    $this->{LineRead}=$LineRead;
    $this->{LineRead}=1;

    $this->{MODELNUM}=1;
    $this->{REPRESENTATIVE_STRUCTURE}=1;

    $this->{LineRead}=0;
    while ($this->{LineRead}<@FileLines){
     if ($this->{FILELINES}[$this->{LineRead}]=~/^HEADER/){
      $this->ReadHeader;}
     elsif ($this->{FILELINES}[$this->{LineRead}]=~/^ATOM/){
      $this->ReadAtom;}
     elsif ($this->{FILELINES}[$this->{LineRead}]=~/^MODEL/){
      $this->ReadModel;}
     $this->{LineRead}++;
    }

    $this->{MODELS}=$this->{MODELNUM};

}#_sub ParseData
#----------------------------------------------------------------------


sub ReadHeader{
    my $this=shift;
    $this->{HEADER}=$this->{FILELINES}[$this->{LineRead}];
    $this->{PDB_ID}=substr($this->{HEADER},62,4);  
    $this->{PDB_CLASSIFICATION}=substr($this->{HEADER},10,30);
    $this->{LAST_MOD_DATE}=substr($this->{HEADER},50,9); 
    $this->{ReadLine}++;
}#_sub ReadHeader

#----------------------------------------------------------------------

sub ReadAtom{

    my $this=shift;

    my $atomline;
    my $resname;  my $resnum; my $resID;  my $resNN;
    my $atomname; my $coordx; my $coordy; my $coordz; my $atomID;
    my $chainID;  my $CHAINID;
    my $altconf;


    $atomline=$this->{FILELINES}[$this->{LineRead}];
    $resname=substr($atomline,17,3);
    $resnum=substr($atomline,22,5);$resnum=~s/ //g;
    $altconf=substr($atomline,16,1);

    $chainID=substr($atomline,21,1);
    if($chainID eq " "){$chainID="a"}; #'a' (no caps!) for atoms without chain ID

    $atomname=substr($atomline,12,4);$atomname=~s/ //g;
    $atomID=$this->{MODELNUM}."_".$chainID."_".$resnum."_".$atomname;
    $resID=$this->{MODELNUM}."_".$chainID."_".$resnum;
    if ($altconf ne " "){$this->{ALTERNATE_CONFORMERS}{$resID}=1}
    $resNN=$resname."_".$resnum;

    $coordx=substr($atomline,30,8);$coordx=~s/ //g;
    $coordy=substr($atomline,38,8);$coordy=~s/ //g;
    $coordz=substr($atomline,46,8);$coordz=~s/ //g;

    $CHAINID="CHAIN_".$chainID;

    if (( not defined $this->{$CHAINID} ) or ($this->{$CHAINID} == 0)){ #[SZB]
	$this->{CHAINNUM}++;
	$this->{CHAINID}[$this->{CHAINNUM}]=$chainID;
	$this->{$CHAINID}=1;
	$this->{CHAIN_N}{$chainID}=$this->{CHAINNUM};
    }

    $this->{FOUNDAMS}{$resname}++;

    # If residue not identified yet, name it and add to sequence
    # and some variables designed to make sequential search easy

     if ($this->{$resID} =~ /^$/){
	$this->{$resID}=$resname;
     	$this->{RESIDUE}[$this->{MODELNUM}][$this->{CHAINNUM}]{$resnum}=$resname;
	$this->{RESNUMS}[$this->{MODELNUM}][$this->{CHAIN_N}{$chainID}].=$resnum."&";
    }    

    $this->{ATOMS}{$resID}.=$atomID."&";
    $this->{FULL_ATOM_LIST}.=$atomID."&";
    $this->{$atomID}[0]=$coordx;
    $this->{$atomID}[1]=$coordy;
    $this->{$atomID}[2]=$coordz;
    $this->{FOUNDATOM}{$atomID}=1;
    $this->{COORDLINE}{$atomID}=$atomline;
    $this->{COORDLINES}{$resID}.="$atomline\n";

}#_sub readAtom



#----------------------------------------------------------------------


sub ReadModel{

    my $this=shift;
    my $modelline=$this->{FILELINES}[$this->{LineRead}];
    my $junk; my $m;
    ($junk,$m)=split(/ +/,$modelline);
    $this->{MODELNUM}=$m;

}#_sub ReadModel


#----------------------------------------------------------------------

# Routine to write structure
# expects: ouptut file (or "STDOUT"), annotation key (see below) model number, 
# chains, residues (see below)

sub WriteStruct{

    my $this=shift;
    my $outfile=shift; 

    if ($outfile ne "STDOUT"){
	open (OF,">$outfile") || die "Cannot create $outfile\n";
    }
    my $annotationkey=shift;

    my $checkmodel=0;
    my $checkres=0;
    my $checkchain=0;

    my $model=shift; if ($model=~/^r/i){$model=$this->{REPRESENTATIVE_STRUCTURE};}
    if ($model ne ""){$checkmodel=1}
    my $chains=shift;  if ($chains ne ""){$checkchain=1}
    my $reslist=shift; if ($reslist ne ""){$checkres=1}

    if (($this->{STRUCTURE_TRANSFORMED}) && (!$this->{COORDS_OK})){$this->V3d2atoms()}

    my $m; my $c; my $cid;
    my $resnums; my @resnums; my $rn; my $resid;
    my $atomids; my @atomids;
    my $ai;
    my $atomline;

    for ($m=1; $m <= $this->{MODELNUM}; $m++){
     if (($checkmodel) && ($m != $model)){next;}
     if ($outfile ne "STDOUT"){printf OF ("MODEL     %4d\n",$m)}
     else {printf ("MODEL     %4d\n",$m)}
     for ($c=1; $c <= $this->{CHAINNUM}; $c++){
	$cid=$this->{CHAINID}[$c];
	if (($checkchain) && ($chains !~ /$cid/)){next;}
	$resnums=$this->{RESNUMS}[$m][$c];$resnums=~s/\&$//;
	@resnums=split(/\&/,$resnums);
	foreach $rn (@resnums){
	    $resid=$m."_".$cid."_".$rn;
	    $atomids=$this->{ATOMS}{$resid};$atomids=~s/\&$//;
	    @atomids=split(/\&/,$atomids);
	    foreach $ai (@atomids){
		$atomline=$this->{COORDLINE}{$ai};
		substr($atomline,30,8)=sprintf("%8.3f",$this->{$ai}[0]);
		substr($atomline,38,8)=sprintf("%8.3f",$this->{$ai}[1]);
		substr($atomline,46,8)=sprintf("%8.3f",$this->{$ai}[2]);
		if ($outfile ne "STDOUT"){print OF "$atomline\n";}
		else {print "$atomline\n";}
	    }#_foreach $ai
	}#_foreach $rn
     }#_for $c
     if ($outfile ne "STDOUT"){print OF "ENDMDL\n";}
     else {print "ENDMDL\n";}

    }#_for $m

    if ($outfile ne "STDOUT"){print OF "END\n";close(OF)}
    else {print "END\n";}

}#_sub WriteStruct

#----------------------------------------------------------------------
# Represent atomic coordinates as a Vector3D objetct
sub AtomVector{

    my $this=shift;
    my $atom1=shift;
    return new Vector3D($this->{$atom1}[0],$this->{$atom1}[1],$this->{$atom1}[2]);

}#sub AtomVector


#----------------------------------------------------------------------

# Getting the vector from Atom1 to Atom2

sub AtomAtomVector{

    my $this=shift;
    my $atom1=shift;
    my $atom2=shift;

    if ($this->{ATOM3D}){
	return Vector3D->Diff2($this->{ATOM_V3D}{$atom1},$this->{ATOM_V3D}{$atom2});
    }
    else {
	return Vector3D->Diff2($this->AtomVector($atom1),$this->AtomVector($atom2));
    }

}#sub AtomAtomVector

#----------------------------------------------------------------------

sub Atoms2v3d{

    my $this=shift;
    my $ai;

    my $fullatomlist=$this->{FULL_ATOM_LIST};
    $fullatomlist=~s/\&$//;
    my @fullatomlist=split(/\&/,$fullatomlist);
    foreach $ai (@fullatomlist){
	$this->{ATOM_V3D}{$ai}=$this->AtomVector($ai);
    }
    $this->{ATOM3D}=1;

}#_sub Atoms2v3d

#----------------------------------------------------------------------

# Routine for test reasons mainly

sub PrintInfo{

    my $this=shift;

    my $c;

    print "$this->{HEADER}\n";
    print "No. of chains found: $this->{CHAINNUM}\n";


}#_sub printInfo

##############################################x


print STDERR"
########################################################
# pdbfile2bmrbnomenclature.pl  < in_pdb  > out_pdb     #
#                                                      #
# This script coverts a pdb file to BMRB nomenclature. #
#                                                      #
# Two-step approach:                                   #
# 1) converting all names to those formally compatible #
#    with BMRB nomenclature (e.g. HB2/3)               #
# 2) assigning correct stereochemistry to geminal      #
#    protons (using an approach based on geometry)     #
#                                                      #
# Will NOT alter bond lengths, angles etc., only       #
# atom NAMES.                                          #
#                                                      #
# Contact: gaspari.zoltan@itk.ppke.hu                  #
########################################################

";

$|=1;
#&Init;
use Getopt::Std;
#use PDBStruct;
#use Vector3D;


getopts('hv') || die "I think not\n";

exit if $opt_h;

$opt_w=lc($opt_w);



# reading PDB file & performing firts conversion step 'on-the-fly'
while(<>){
    chomp;
    # converting everything what we can here...
    if ($_ =~/^ATOM|^HETATM/){
	$resname=substr($_,17,3);
        # As a rule, will convert all geminal protons with label 1 to 3
        # stereospecific labels will be adjusted at stage 2 anyway
        # and in this way we do not rely on a dictionary but handle everything (hopefully)
        # NOT YET COMPLETE
	if ($resname =~/ARG|ASP|ASN|CYS|GLU|GLN|HIS|LEU|LYS|MET|PHE|PRO|SER|TRP|TYR/){
	    $_=~s/ [13]HB / HB3 /;
	    $_=~s/ HB1 / HB3 /;
	    $_=~s/ 2HB / HB2 /;
	}
	if ($resname =~/ARG|GLU|GLN|LYS|MET|PRO/){
	    $_=~s/ [13]HG / HG3 /;
	    $_=~s/ HG1 / HG3 /;
	    $_=~s/ 2HG / HG2 /;
	}
	if ($resname =~/ARG|LYS|MET|PRO/){
	    $_=~s/ [13]HD / HD3 /;
	    $_=~s/ HD1 / HD3 /;
	    $_=~s/ 2HD / HD2 /;
	}
	if ($resname =~/LYS/){
	    $_=~s/ [13]HE / HE3 /;
	    $_=~s/ HE1 / HE3 /;
	    $_=~s/ 2HE / HE2 /;
	    $_=~s/ ([123])HZ / HZ$1 /;$_=~s/ HNZ([123]) / HZ$1  /; # check the spaces..
	}
	if ($resname =~/ASN/){
	    $_=~s/ ([12])HD2 / HD2$1 /;
	}
	if ($resname =~/GLN/){
	    $_=~s/ ([12])HE2 / HE2$1 /;
	}
	if ($resname =~/THR/){
	    $_=~s/ ([123])HG2 / HG2$1 /;
	}
	if ($resname =~/VAL/){
	    $_=~s/ ([123])HG1 / HG1$1 /;
	    $_=~s/ ([123])HG2 / HG2$1 /;
	}
	if ($resname =~/LEU/){
	    $_=~s/ ([123])HD1 / HD1$1 /;
	    $_=~s/ ([123])HD2 / HD2$1 /;
	}
	if ($resname =~/ILE/){
	    $_=~s/ CD  / CD1 /;
	    $_=~s/ ([123])HG2 / HG2$1 /;
	    $_=~s/ ([13])HG1 / HG13 /;
	    $_=~s/ HG1([13]) / HG13 /;
	    $_=~s/ 2HG1 / HG12 /;
	    $_=~s/ ([123])HD1 / HD1$1 /;
	    $_=~s/ HD1([123]) / HD1$1 /;
	    $_=~s/  ([123])HD / HD1$1 /;
	    $_=~s/  HD([123]) / HD1$1 /;
	}
	if ($resname =~/GLY/){
	    $_=~s/ [13]HA / HA3 /;
	    $_=~s/ HA1 / HA3 /;
	    $_=~s/ 2HA / HA2 /;
	}
        # This also seems necessary...
	if ($resname =~/ARG/){
	    $_=~s/ 1HH1 / HH11 /;
	    $_=~s/ 2HH1 / HH12 /;
	    $_=~s/ 1HH2 / HH21 /;
	    $_=~s/ 2HH2 / HH22 /;
	}
  
        # Amide HN to H - only for amino acids!
	$_=~s/ HN / H  / unless ($_=~/^HETATM/);

    }

    # passing the lines to a PDStruct object for the 2nd stage    
    $structlines.="$_\n";
}
$in_pdb=new PDBStruct($structlines);

# code from PDBStruct.pm (routine WriteStruct)

for ($m=1; $m <= $in_pdb->{MODELNUM}; $m++){
    #print "MODEL: $m chainnum: $in_pdb->{MODELNUM}\n";
     #if (($checkmodel) && ($m != $model)){next;}
     for ($c=1; $c <= $in_pdb->{CHAINNUM}; $c++){
        #print "CHAIN: $c\n";
        $cid=$in_pdb->{CHAINID}[$c];
        if (($checkchain) && ($chains !~ /$cid/)){next;}
        $resnums=$in_pdb->{RESNUMS}[$m][$c];$resnums=~s/\&$//;
        @resnums=split(/\&/,$resnums);
        foreach $rn (@resnums){
            $resid=$m."_".$cid."_".$rn;
	     if ($in_pdb->{$resid} eq "ARG"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","CD","CG","HG2","HG3");
		&checkgeminal($resid,"CG","NE","CD","HD2","HD3");
	     }
	     if ($in_pdb->{$resid} eq "ASP"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "ASN"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "CYS"){
		&checkgeminal($resid,"CA","SG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "GLU"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","CD","CG","HG2","HG3");
	     }
	     if ($in_pdb->{$resid} eq "GLN"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","CD","CG","HG2","HG3");
	     }
	     if ($in_pdb->{$resid} eq "GLY"){
		&checkgeminal($resid,"N","C","CA","HA2","HA3");
	     }
	     if ($in_pdb->{$resid} eq "HIS"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "ILE"){
		&checkgeminal($resid,"CB","CD1","CG1","HG12","HG13");
	     }
	     if ($in_pdb->{$resid} eq "LEU"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "LYS"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","CD","CG","HG2","HG3");
		&checkgeminal($resid,"CG","CE","CD","HD2","HD3");
		&checkgeminal($resid,"CD","NZ","CE","HE2","HE3");
	     }
	     if ($in_pdb->{$resid} eq "MET"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","SD","CG","HG2","HG3");
	     }
	     if ($in_pdb->{$resid} eq "PHE"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "PRO"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
		&checkgeminal($resid,"CB","CD","CG","HG2","HG3");
		&checkgeminal($resid,"CG","N" ,"CD","HD2","HD3");
	     }
	     if ($in_pdb->{$resid} eq "SER"){
		&checkgeminal($resid,"CA","OG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "TRP"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }
	     if ($in_pdb->{$resid} eq "TYR"){
		&checkgeminal($resid,"CA","CG","CB","HB2","HB3");
	     }

	      if ($in_pdb->{$resid} eq "VAL"){
		$t=checkgeminal($resid,"CA","HB","CB","CG1","CG2");
		if ($t < 0){
		    for ($x=1; $x<=3; $x++){
			$Xid=$resid."_HG1".$x;
			$Yid=$resid."_HG2".$x;
			$tmpx=$in_pdb->{$Xid}[0];
			$tmpy=$in_pdb->{$Xid}[1];
			$tmpz=$in_pdb->{$Xid}[2];
			$in_pdb->{$Xid}[0]=$in_pdb->{$Yid}[0];
			$in_pdb->{$Xid}[1]=$in_pdb->{$Yid}[1];
			$in_pdb->{$Xid}[2]=$in_pdb->{$Yid}[2];
			$in_pdb->{$Yid}[0]=$tmpx;
			$in_pdb->{$Yid}[1]=$tmpy;
			$in_pdb->{$Yid}[2]=$tmpz;
		    }
		    print STDERR "Methyl swap performed for $resid\n" if $opt_v;
		}
	      }
	      if ($in_pdb->{$resid} eq "LEU"){
		$t=checkgeminal($resid,"CA","HB","CB","CD1","CD2");
		if ($t < 0){
		    for ($x=1; $x<=3; $x++){
			$Xid=$resid."_HD1".$x;
			$Yid=$resid."_HD2".$x;
			$tmpx=$in_pdb->{$Xid}[0];
			$tmpy=$in_pdb->{$Xid}[1];
			$tmpz=$in_pdb->{$Xid}[2];
			$in_pdb->{$Xid}[0]=$in_pdb->{$Yid}[0];
			$in_pdb->{$Xid}[1]=$in_pdb->{$Yid}[1];
			$in_pdb->{$Xid}[2]=$in_pdb->{$Yid}[2];
			$in_pdb->{$Yid}[0]=$tmpx;
			$in_pdb->{$Yid}[1]=$tmpy;
			$in_pdb->{$Yid}[2]=$tmpz;
		    }
		    print STDERR "Methyl swap performed for $resid\n" if $opt_v;
		}
	      }

	}
     }
}

$in_pdb->WriteStruct("STDOUT");

exit;

sub checkgeminal{
    my $rid=shift;
    my $down=shift; # atom closer to bb
    my $up=shift; # atom further from bb
    my $front=shift; # atom with the geminal Hs
    my $X=shift; # geminal H 
    my $Y=shift; # geminal H
    $downid=$rid."_".$down;
    $upid=$rid."_".$up;
    $frontid=$rid."_".$front;
    $Xid=$rid."_".$X;
    $Yid=$rid."_".$Y;

    # principle: compare the vectorial product 
    # of the X-Y and vectors with the 
    # sum of X-front and Y-front vectors to the direction 
    # of the down-up vector, and go from there.

    $du_vector=$in_pdb->AtomAtomVector($downid,$upid);
    $XY_vector=$in_pdb->AtomAtomVector($Xid,$Yid);
    $front_vector=$in_pdb->AtomAtomVector($Xid,$frontid)->Sum($in_pdb->AtomAtomVector($Yid,$frontid));
    $frontXY_vp=$front_vector->VectorialProduct($XY_vector);
    # Normalizing...
    $frontXY_vp_norm=$frontXY_vp->Norm();
    $du_vector_norm=$du_vector->Norm();
    $test=$frontXY_vp_norm->ScalarProduct($du_vector_norm);
    printf (">> $resid | $X <-> $Y | Test: %7.3f\n",$test) if $opt_v;

    if ($test < 0){
	$tmpx=$in_pdb->{$Xid}[0];
	$tmpy=$in_pdb->{$Xid}[1];
	$tmpz=$in_pdb->{$Xid}[2];
	$in_pdb->{$Xid}[0]=$in_pdb->{$Yid}[0];
	$in_pdb->{$Xid}[1]=$in_pdb->{$Yid}[1];
	$in_pdb->{$Xid}[2]=$in_pdb->{$Yid}[2];
	$in_pdb->{$Yid}[0]=$tmpx;
	$in_pdb->{$Yid}[1]=$tmpy;
	$in_pdb->{$Yid}[2]=$tmpz;
    }

    return $test;
}


























