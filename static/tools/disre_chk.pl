#!/usr/bin/perl 

#To elimninate dependence on an external package, Vector3D is defined here
{

package Vector3D;

###################################################
##
## Class for handling 3D points and vectors
## The constructor expects an array with the coordinates of
## the point. 
##
## All subroutines are also defined in a static manner,
## e.g Length() and Length2() differ in invocation
## (the first is invoked as an object method, the second
## as a class method with an object as argument)
##
## Parts of the code are based on the following book:
## [JT]: L. Jánossy, P. Tasnádi: Vektorszámítás (1982)
## Tankönyvkiadó, Budapest (in Hungarian)
##
## Version 0.2, 13.02.2002. 
###################################################

use strict;


##############
## Constructor
##############

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

#    $this->{LENGTH}=sqrt(($this->{COORDX}**2)+($this->{COORDY}**2)+($this->{COORDZ}**2));
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

# Making unity vector from $this

sub Norm2{

    my $type1=shift;
    my $first=shift;
    return ($first->Norm());

}#_sub Norm2


#-----------------------------------------------------

# Scaling $this (with argument $scale)

sub Scale{

    my $this=shift;
    my $scale=shift;
    my $length=$this->Length();
    my $m_x=$this->{COORDX}*$scale;
    my $m_y=$this->{COORDY}*$scale;
    my $m_z=$this->{COORDZ}*$scale;

    return (new Vector3D($m_x,$m_y,$m_z));

}#_sub Scale

#-----------------------------------------------------

# Scaling Vector3D object (fisrt arg) with $scale (second arg)

sub Scale2{

    my $type1=shift;
    my $first=shift;
    my $scale=shift;
    return ($first->Scale($scale));

}#_sub Scale2

#-----------------------------------------------------

# Distance calculator for $this and another Vector3D object

sub Distance{

    my $this=shift;
    my $second=shift;

    my $xd=$this->{COORDX}-$second->{COORDX};
    my $yd=$this->{COORDY}-$second->{COORDY};
    my $zd=$this->{COORDZ}-$second->{COORDZ};

    my $dist=sqrt(($xd**2)+($yd**2)+($zd**2));

    return($dist);
}#_sub Distance

#-----------------------------------------------------

# Distance calculator for two Vector3D objects

sub Distance2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;

    return($first->Distance($second));

}#_sub Distance2

#-------------------------------------------

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


#--------------------------------------------

# Harmas vegyes szorzat : a*(b X c)

sub MixedProduct{

    my $this=shift;
    my $second=shift;
    my $third=shift;

    return ($this->ScalarProduct($second->VectorialProduct($third)));

}#_sub MixedProduct


#-------------------------------------------

# Harmas vegyes szorzat 2

sub MixedProduct2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;
    my $third=shift;

    return ($first->MixedProduct($second,$third));


}#_sub MixedProduct2

#--------------------------------------------
sub MatrixMul{

   my $this=shift;
   my $row1=shift;
   my $row2=shift;
   my $row3=shift;

   my $x1=$this->ScalarProduct($row1);
   my $x2=$this->ScalarProduct($row2);
   my $x3=$this->ScalarProduct($row3);

   return new Vector3D($x1,$x2,$x3);


}#_sub MatrixMul

#--------------------------------------------

# Reciprocal vector of $this (needs two arguments)
# [JT, p. 4., eq. 4.8-10.]

sub Reciprocal{

    my $this=shift;
    my $second=shift;
    my $third=shift;

    my $v=$this->MixedProduct($second,$third);
#   print "[Reciprocal] $v";

    my $rec=$second->VectorialProduct($third);
    $rec=$rec->Scale((1/$v));
#   $rec->PrintInfo();
    return ($rec);


}#_sub MixedProduct


#-------------------------------------------

# Reciprocal vector, three args

sub Reciprocal2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;
    my $third=shift;

    return ($first->Reciprocal($second,$third));

}#_sub Reciprocal2

#--------------------------------------------

# Projecting $this to a plane defined by its normal vector
# A point in the plane ($a) and its normal vector ($norm) 
# is prvided as an argument
# [JT, pp. 59] $M=$normal; $b=$this

sub Project{

    my $this=shift;
    my $a=shift;
    my $normal=shift;

    # The projaction line is perpendicular to the plane, so $M=$normal
    # $this is the point to project, it is on the line
    # Ensuring consistency with symbols in [JT]
    my $M=$normal;  
    my $b=$this;

    # Constructing two vectors in the plane:
    # A vector not paralel to $norm
    my $n2=new Vector3D($normal->{COORDY},((-1)*$normal->{COORDX})+5,$normal->{COORDZ});
    # $K is perpendicular to $normal, $L is perpendicular to $normal and $K
    my $K=$normal->VectorialProduct($n2);
    my $L=$normal->VectorialProduct($K);

#    print "K,L,M:\n";
#    $K->PrintInfo();
#    $L->PrintInfo();
#    $M->PrintInfo();
#    my $mixed=$K->MixedProduct($L,$M);
#    print "$mixed\n\n";
#    my $vecKL=$K->VectorialProduct($L);
#    $vecKL->PrintInfo();
#    my $svecKL=$vecKL->Scale((1/$mixed));
#    $svecKL->PrintInfo();


    # Constructing reciprocal vectors
    my $k=$K->Reciprocal($L,$M);
    my $l=$L->Reciprocal($M,$K);
    my $m=$M->Reciprocal($K,$L);

#    print "k,l,m, vegyesszorzatuk:\n";
#    $k->PrintInfo();
#    $l->PrintInfo();
#    $m->PrintInfo();
#    $mixed=$k->MixedProduct($l,$m);
#    print "$mixed\n\n";

    my $s=$m->ScalarProduct($b->Diff($a));

    # $r is the intercept of the line and the plane,
    # in other words, the projected point.
    my $r=$b->Diff($M->Scale($s));

    return($r);

}#_sub Project

#--------------------------------------------

sub Project2{

    my $type1=shift;
    my $b=shift;
    my $a=shift;
    my $normal=shift;

    return($b->Project($a,$normal));

}#_sub Project2

#---------------------------------------------

# Angle of $this and a second Vector3D object

sub Angle{

    my $this=shift;
    my $second=shift;

    my $angle;

    my $this_length=$this->Length();
    my $second_length=$second->Length();

    my $scalar=$this->ScalarProduct($second);
    my $vectorial=$this->VectorialProduct($second);

    if ($this_length*$second_length == 0){$angle="NaN"}
    else {

	my $cos=$scalar/($this_length*$second_length);
	my $sin=$vectorial->Length()/($this_length*$second_length);

#       print STDERR "$scalar $vectorial $this_length $second_length $sin $cos\n";

	$angle=atan2($sin,$cos);
    }

    return($angle);

}#_sub Angle


#-------------------------------------------

# Angle of two Vector3D objects

sub Angle2{

    my $type1=shift;
    my $first=shift;
    my $second=shift;

    return ($first->Angle($second));

}#_sub Angle2

#--------------------------------------------

# Rotation around X axis
# Expects angle (radian)
# Currently not uses rotation matrices (mainly for
# "speed" reasons :-))

sub RotateX{

    my $this=shift;
    my $angle=shift;
    my $sin=sin($angle);
    my $cos=cos($angle);

    my $resultx; my $resulty; my $resultz;

    $resultx=$this->{COORDX};
    $resulty=$cos*$this->{COORDY}-$sin*$this->{COORDZ};
    $resultz=$sin*$this->{COORDY}+$cos*$this->{COORDZ};

    return new Vector3D($resultx,$resulty,$resultz);

}#_sub RotateX

#-------------------------------------------

# Rotation around Y axis
# Expects angle (radian)
# Currently not uses rotation matrices (mainly for
# "speed" reasons :-))

sub RotateY{

    my $this=shift;
    my $angle=shift;
    my $sin=sin($angle);
    my $cos=cos($angle);

    my $resultx; my $resulty; my $resultz;

    $resultx=$cos*$this->{COORDX}+$sin*$this->{COORDZ};
    $resulty=$this->{COORDY};
    $resultz=$cos*$this->{COORDZ}-$sin*$this->{COORDX};

    return new Vector3D($resultx,$resulty,$resultz);

}#_sub RotateY

#-------------------------------------------

# Rotation around Z axis
# Expects angle (radian)
# Currently not uses rotation matrices (mainly for
# "speed" reasons :-))

sub RotateZ{

    my $this=shift;
    my $angle=shift;
    my $sin=sin($angle);
    my $cos=cos($angle);

    my $resultx; my $resulty; my $resultz;

    $resultx=$cos*$this->{COORDX}-$sin*$this->{COORDY};
    $resulty=$sin*$this->{COORDX}+$cos*$this->{COORDY};
    $resultz=$this->{COORDZ};

    return new Vector3D($resultx,$resulty,$resultz);

}#_sub RotateZ

#-------------------------------------------

# Rotation around X, Y and Z axis (in this order)
# Expects 3 angles (radian)
# Currently not uses rotation matrices (mainly for
# "speed" reasons :-))

sub Rotate{

    my $this=shift;
    my $anglex=shift;
    my $angley=shift;
    my $anglez=shift;

    # Mmmm.... piece of cake :-)))
    return (($this->RotateX($anglex))->RotateY($angley))->RotateZ($anglez);
#    return $this->RotateZ($this->RotateY($this->RotateX($anglex),$angley),$anglez);

}#_sub Rotate

#-------------------------------------------


}


$|=1;

print STDERR "
                     :-] disre_chk.pl [-:
                 with a GROMACS-like interface
                   version as of 27 Jan 2017

Option     Filename  Type          Description
------------------------------------------------------------
";


use Getopt::Std;
#use Vector3D;

$opt_f="eiwit.pdb";
$opt_d="noe.disre";
$opt_o="disre_chk.out";
$opt_U="unviolated.disre";
$opt_V="violated.disre";
$opt_a=-6;
$opt_t=0;
$opt_e=-6;

getopts('f:d:o:U:V:a:e:t:rvh') || die "Please check your options\n";

if ($opt_h){$HO="yes"} else{$HO="no"}
if ($opt_r){$RO="yes"} else{$RO="no"}
if($opt_v){$V=1;$VO="yes"}else{$VO="no"}

printf STDERR ("  -f  %12s   Input         generic structure file (multi-model PDB at the moment\n",$opt_f);
printf STDERR ("  -d  %12s   Input         disre file\n",$opt_d);
printf STDERR ("  -o  %12s   Output        NOE check report file\n",$opt_o);
printf STDERR ("  -U  %12s   Output        disre file with unviolated restraints\n",$opt_U);
printf STDERR ("  -V  %12s   Output        disre file with violated restraints\n",$opt_V);

print STDERR"
      Option   Type  Value  Description
------------------------------------------------------
";

printf STDERR ("          -h   bool %9s  Print help info and quit\n",$HO);
printf STDERR ("          -v   bool %9s  Be loud and noisy\n",$VO);
printf STDERR ("          -a   int  %9s  Exponent for ambiguous restraint averaging\n",$opt_a);
printf STDERR ("          -e   int  %9s  Exponent for ensemble averaging\n",$opt_e);
printf STDERR ("          -t   float  %5.2f  Tolerance for violation\n",$opt_t);
printf STDERR ("          -r   bool %9s  Rename hydrogens like 2HB -> HB2 in PDB file\n\n",$RO);


if ($opt_h) {
die "

Trying to do meaningful NOE analysis on a bunch of structures. 
No warranty.

";
}

$tolerance=$opt_t;
open (TF,"$opt_d") || die "Cannot open disre file $opt_d for reading\n";
$readdisre=0;
$disrenum=-1;
$pairnum=0;
$previndex=-1;
$pdbok=1;
$violated=0;
while (<TF>){
    chomp;
    next if ($_=~/^ *\#/);
    if ($_=~/^ *[0-9]+/){

       ($disre,$comment)=split(/\#/,$_);
       #if ($comment !~ /[A-Z]/){$pdbok=0}
       $disre=~s/^ +//;
       ($index,$rnum1,$rname1,$atomname1,$rnum2,$rname2,$atomname2,$up1)=split(/ +/,$disre);
       #($ai,$aj,$type,$index,$typea,$low,$up1,$up2,$fac)=split(/ +/,$disre);
       #($atom1,$atom2)=split(/\:/,$comment);
       #$atom1=~s/^ *([0-9]+) *([A-Z][A-Z0-9]*) *$/$1\_$2/;
       #$atom2=~s/^ *([0-9]+) *([A-Z][A-Z0-9]*) *$/$1\_$2/;
       #print "$_ :: $comment :: $atom1 <> $atom2\n" if $V;
       $atom1=$rnum1."_".$atomname1;
       $atom2=$rnum2."_".$atomname2;

       # Storing atoms to read
       $atomtoread{$atom1}=1;
       $atomtoread{$atom2}=1;

       $PAIR[$pairnum]=$atom1."|".$atom2;

       # Keeping intrenal disre numbering with a reference to the original one
       if ($index != $previndex){$disrenum++;$INDEX[$disrenum]=$index} 
       $PAIRSFORDR[$disrenum].=":".$pairnum;

       # Checking & bookkeeping
       if (($RESNAME{$rnum1} =~ /[A-Z]/) && ($RESNAME{$rnum1} ne $rname1)){
 	die "In DR $dr residue $rnum1 is a $rname1 but was prevously identified a $RESNAME{$rnum1}\n";
       }
       else{
	$RESNAME{$rnum1}=$rname1;
       }
       if (($RESNAME{$rnum2} =~ /[A-Z]/) && ($RESNAME{$rnum2} ne $rname2)){
	die "In DR $dr residue $rnum2 is a $rname2 but was prevously identified a $RESNAME{$rnum2}\n";
       }
       else{
	$RESNAME{$rnum2}=$rname2;
       }

       $DIST[$disrenum]=$up1;

       print "$_ | $index DR $disrenum PR $pairnum PFDR $PAIRSFORDR[$disrenum]\n" if $V;
       $pairnum++;
       $previndex=$index;
    }

}#_while TF
   
close(TF);
# numbering starts with 0...
$dno=$disrenum+1;
print STDERR "Read $pairnum atom pairs in $dno distance restraints\n";
if (!$pdbok){
    print STDERR "Warning: one or more restraint lines contain invalid atom names, treating PDB file as input may fail\n";
}

# Guessing input structure file format from extension
if ($opt_f =~ /\.pdb$|\.ent$/){$SFORMAT="pdb"}
else {
    die "Unrecognized struture file format, currently only PDB files (.pdb, .ent) are accepted\n";
}


open (PDB,"$opt_f") || die "Cannot open PDB file $opt_f\n";
$MODEL=0; 
while (<PDB>){
  chop;
  if ($_ =~ /^MODEL/){$MODEL++}
  elsif ($_=~/^ATOM/) {

   $resnum=substr($_,22,5);$resnum=~s/ //g;
   $name=substr($_,12,4);$name=~s/ //g;

   #if ($name eq "H"){$name="HN"}
   $name=~s/([0-9])H([A-Z0-9]+)/H$2$1/ if $opt_r;

   $atomexp=$resnum."_".$name;
   #print "$_ :: $atomexp\n";
 
   if ($atomtoread{$atomexp}) {
       $coordx=substr($_,30,8);$coordx=~s/ //g;
       $coordy=substr($_,38,8);$coordy=~s/ //g;
       $coordz=substr($_,46,8);$coordz=~s/ //g;
       $atomv{$MODEL}{$atomexp}=new Vector3D($coordx,$coordy,$coordz);
       $foundatom{$atomexp}=1;
       #print "  \#MODEL $MODEL atomexp $atomexp\#";$atomv{$MODEL}{$atomexp}->PrintInfo();
   }#_if
 }#_if ATOM
}#_while PDB
close(PDB);
# And now, more music

open (OF,">$opt_o") || die "Cannot create output file $opt_o\n";
open (UF,">$opt_U") || die "Cannot create output file $opt_U\n";
open (VF,">$opt_V") || die "Cannot create output file $opt_V\n";
$unvioln=0;
for ($drn=0; $drn <= $disrenum; $drn++){
    $pairs=$PAIRSFORDR[$drn];
    $pairs=~s/^://;
    @pr=split(/\:/,$pairs);
    $pn=@pr;

    $ensembleavedist=0;
    @AAVE=();

    for ($M=1; $M <= $MODEL; $M++){
	$ambigavedist=0;
	printf OF ("# DR %4d MODEL %3d",$drn,$M) if ($pn > 1);
	foreach $pair (@pr){
	($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
	printf OF (" | PAIR $pair $atomid1 $atomid2") if ($pn > 1);
            if ($M == 1){
	     #print STDERR "  \#MODEL $M atom1 $atomid1\#";$atomv{$M}{$atomid1}->PrintInfo() if $opt_v;
	     #print STDERR "  \#MODEL $M atom2 $atomid2\#";$atomv{$M}{$atomid2}->PrintInfo() if $opt_v;
            }
            if ($foundatom{$atomid1}*$foundatom{$atomid2}){
	     $dist=($atomv{$M}{$atomid1}->Diff($atomv{$M}{$atomid2}))->Length();
	     $ambigavedist+=exp($opt_a*log($dist));
             printf OF (" %5.3f",$dist) if ($pn > 1); 
           }
           if (!$foundatom{$atomid1}){print STDERR "ATOM $atomid1 not found\n"}
           if (!$foundatom{$atomid2}){print STDERR "ATOM $atomid2 not found\n"}
	}
	$ambigavedist=exp((1/$opt_a)*log($ambigavedist/$pn));
	$AAVE[$M]=$ambigavedist;
        printf OF (" aaved=%5.3f\n",$ambigavedist) if ($pn > 1); 
	$ensembleavedist+=exp($opt_e*log($ambigavedist));
    }
    $ensembleavedist=exp((1/$opt_e)*log($ensembleavedist/$MODEL));
    printf OF ("%5d %5.2f  <> ",$drn,$DIST[$drn]);
    for ($M=1; $M <= $MODEL; $M++){
        printf OF (" %5.3f",$AAVE[$M]);
    }
    #printf OF ("   %5.3f\n",$ensembleavedist);

    $violation=0;
    if ($ensembleavedist > ($DIST[$drn]+$tolerance)){
	$violated++;
	#$violation=$ensembleavedist-($DIST[$drn]+$tolerance);
	$violation=$ensembleavedist-($DIST[$drn]);
        #print STDERR "violation: $violation\n";
	#$aveviol+=$ensembleavedist-$DIST[$drn];
	$aveviol+=$violation;
	#if (($ensembleavedist-$DIST[$drn]) > $maxviol){
	#    $maxviol=$ensembleavedist-$DIST[$drn];
	#    $maxvioldr=$drn;
	#}
	if ($violation > $maxviol){
	    $maxviol=$violation;
	    $maxvioldr=$drn;
	}
        # writing violated list
	printf VF ("#Was $INDEX[$drn] in original list | ACTUAL %7.3f \n",$ensembleavedist);
	foreach $pair (@pr){
	  ($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
          ($rnum1,$atom1)=split(/\_/,$atomid1);
          ($rnum2,$atom2)=split(/\_/,$atomid2);
          # ensuring for output disre list to start with 1 (was violated-1 in the original code)
	  printf VF (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$violated,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$drn]);
	}
    }
    # writing unviolated list
    else{
	printf UF ("#Was $INDEX[$drn] in original list | ACTUAL %7.3f \n",$ensembleavedist);
	foreach $pair (@pr){
	  ($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
          ($rnum1,$atom1)=split(/\_/,$atomid1);
          ($rnum2,$atom2)=split(/\_/,$atomid2);
          # +1 added for output disre list to start with 1
	  printf UF (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$unvioln+1,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$drn]);
	}
	$unvioln++;
    }
    printf OF (" EAD %5.3f VIOL %5.3f\n",$ensembleavedist,$violation);

}
if ($violated > 0){
    $aveviol/=$violated;
}
printf OF ("# There are %d violated restraints. Average violation is %5.3f A, maximum violation is %5.3f A",$violated,$aveviol,$maxviol);
if ($maxviol){print OF " (restraint $INDEX[$maxvioldr])"}
print OF "\n";

printf STDERR ("# There are %d violated restraints. Average violation is %5.3f A, maximum violation is %5.3f A",$violated,$aveviol,$maxviol);
if ($maxviol){print STDERR " (restraint $INDEX[$maxvioldr])"}
print STDERR "\n";

close(OF);
close(UF);
close(VF);
