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

find_aa_bb.pl - script to find sites in Mom and Dad's bam files where Mom is AA and Dad is BB.

=head1 SYNOPSIS

  find_aa_bb.pl --mombam <mother's BAM file> --dadbam <father's BAM file> --ref <reference fasta> --min <minimum coverage required in each BAM>

=head1 DESCRIPTION

This script reports sites for which all reads from Mom are one allele and all reads from Dad are the other, and both Mom and Dad have at least the minimum required number of reads.

=cut

#------------
# Begin MAIN
#------------

process_commandline();

$ENV{'PATH'} = "/home/nhansen/projects/bard/scripts:/home/nhansen/projects/bard/c:/home/nhansen/projects/bamcounts:$ENV{PATH}";

my $mombam = $Opt{'mombam'};
my $dadbam = $Opt{'dadbam'};
my $minreads = $Opt{'min'};
my $ref = $Opt{'ref'};

read_bam_files($mombam, $dadbam, $ref, $minreads);

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( min => 20 , countopts => '' );
    GetOptions(
        \%Opt, qw(
            manual help+ version
            verbose mombam=s dadbam=s ref=s min=i
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "find_aa_bb.pl, ", q$Revision: 4129 $, "\n";
	}

	if ( !$Opt{'mombam'} || !$Opt{'dadbam'}) {
		pod2usage( 0 );
                die; 
	}

	if ( !exists( $Opt{ref} ) ) {
		die "Specify reference fasta filename with option --ref.\n";
	}
}

sub read_bam_files {
    my $mombam = shift;
    my $dadbam = shift;
    my $ref = shift;
    my $minreads = shift;

    open MOM, "bamcounts -fasta $ref -bam $mombam $Opt{countopts} | "
        or die "Couldn\'t run bamcounts on Mom!\n";

    open DAD, "bamcounts -fasta $ref -bam $dadbam $Opt{countopts} | "
        or die "Couldn\'t run bamcounts on Dad!\n";

    my $rh_mom_homs = {};
    my $rh_dad_homs = {};
    my ($mom_line, $dad_line);
    while (defined($mom_line = <MOM>) || defined($dad_line = <DAD>)) {
        if ($mom_line && $mom_line =~ /^(\S+)\s(\d+)\s(\S)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($chr, $pos, $ref, $A_count, $T_count, $G_count, $C_count, $a_count, $t_count, $g_count, $c_count) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
            my %mom_base_count;
            $mom_base_count{'A'} = $A_count + $a_count;
            $mom_base_count{'C'} = $C_count + $c_count;
            $mom_base_count{'G'} = $G_count + $g_count;
            $mom_base_count{'T'} = $T_count + $t_count;
            # how many bases represented?
            my @bases = grep { $mom_base_count{$_}>0 } keys %mom_base_count;
            if (@bases==1) {
                my $mom_base = $bases[0];
                if ($mom_base_count{$mom_base} >= $minreads) { # record
                    if ($rh_dad_homs->{"$chr:$pos"}) {
                        if ($rh_dad_homs->{"$chr:$pos"}->{'base'} ne $mom_base) {
                            my $posmo = $pos - 1;
                            print "$chr\t$posmo\t$pos\n"; 
                        }
                        delete $rh_dad_homs->{"$chr:$pos"};
                    }
                    else { # save for later
                        $rh_mom_homs->{"$chr:$pos"}->{'base'} = $mom_base;
                        #print "MOM: $chr:$pos ($mom_base)\n";
                    }
                }
            }
            $dad_line = <DAD>;
        }
        if ($dad_line && $dad_line =~ /^(\S+)\s(\d+)\s(\S)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($chr, $pos, $ref, $A_count, $T_count, $G_count, $C_count, $a_count, $t_count, $g_count, $c_count) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
            my %dad_base_count;
            $dad_base_count{'A'} = $A_count + $a_count;
            $dad_base_count{'C'} = $C_count + $c_count;
            $dad_base_count{'G'} = $G_count + $g_count;
            $dad_base_count{'T'} = $T_count + $t_count;
            # how many bases represented?
            my @bases = grep { $dad_base_count{$_}>0 } keys %dad_base_count;
            if (@bases==1) {
                my $dad_base = $bases[0];
                if ($dad_base_count{$dad_base} >= $minreads) { # record
                    if ($rh_mom_homs->{"$chr:$pos"}) {
                        if ($rh_mom_homs->{"$chr:$pos"}->{'base'} ne $dad_base) {
                            my $posmo = $pos - 1;
                            print "$chr\t$posmo\t$pos\n"; 
                        }
                        delete $rh_mom_homs->{"$chr:$pos"};
                    }
                    else { # save for later
                        $rh_dad_homs->{"$chr:$pos"}->{'base'} = $dad_base;
                        #print "DAD: $chr:$pos ($dad_base)\n";
                    }
                }
            }
        }
    }
}

close MOM;
close DAD;

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--minscore> N

Minimum number of bases required in each of the two samples.

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
