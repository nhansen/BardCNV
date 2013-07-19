#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr( q$Rev: 0$, 4 ) / 10000;

=head1 NAME

calc_ploidy.pl - script to examine a model file and report ploidy based on optimal values of mu and contam.

=head1 SYNOPSIS

  calc_ploidy.pl <modelfile> <normalbam> <tumorbam>

=head1 DESCRIPTION

This script reads in the file of optimal parameters from BARD, and calculates estimated purity and ploidy from the values of mu and rho_contam.

=cut

#------------
# Begin MAIN
#------------

my $Usage = qq!calc_ploidy.pl <modelfile> <normalbam> <tumorbam> <train error file>\n!;

$#ARGV>=2
    or die "$Usage";

my $modelfile = $ARGV[0];
my $normalbam = $ARGV[1];
my $tumorbam = $ARGV[2];
my $trainerror = $ARGV[3];

open MODEL, $modelfile
    or die "Couldn\'t open $modelfile: $!\n";

my ($mu, $rhocontam);
while (<MODEL>) {
    if (/^mu= (\S+)$/) {
        $mu = $1;
    }
    elsif (/^rho_contam= (\S+)$/) {
        $rhocontam = $1;
    }
}
close MODEL;

if (!$mu) {
    die "Failed to obtain value of mu from $modelfile!\n";
}

if (!$rhocontam) {
    die "Failed to obtain value of rho_contam from $modelfile!\n";
}

my $normal_reads = `samtools idxstats $normalbam | awk '{sum += \$3} END {print sum}'`;
chomp $normal_reads;

if (!$normal_reads) {
    print "Couldn\'t obtain number of normal reads from file $normalbam!\n";
}

my $tumor_reads = `samtools idxstats $tumorbam | awk '{sum += \$3} END {print sum}'`;
chomp $tumor_reads;

if (!$tumor_reads) {
    print "Couldn\'t obtain number of tumor reads from file $tumorbam!\n";
}

my $ploidy = 1/(1-$rhocontam) * ($tumor_reads/($mu*$normal_reads) - $rhocontam*2);

my $purity = 1 - $rhocontam;

my $mll;
if ($trainerror) {
    $mll = `grep -i mll $trainerror | tail -1 | awk '{print \$NF}' | sed 's/)//'`;
    chomp $mll;
}

print "Purity $purity, ploidy $ploidy for MU $mu RHO $rhocontam MLL $mll\n";

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--normal_counts>

Specify the name of an output file from "bamcounts" for the normal sample.

=item B<--tumor_counts>

Specify the name of an output file from "bamcounts" for the tumor sample.

=item B<--mindepth>

Specify the minimum total read depth (total of normal and tumor) to report a site.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use.

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose.

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation.

=cut
