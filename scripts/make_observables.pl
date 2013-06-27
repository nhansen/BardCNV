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

make_observables.pl - script to write out positions with read depth ratio and B-allele frequency values for use in training and CNV state prediction with bardCNV.

=head1 SYNOPSIS

  make_observables.pl --normal_counts <counts from normal file> --tumor_counts <counts from tumor file>

=head1 DESCRIPTION

This script reads in the file of counts from the normal sample, then, based on actual counts or a lack of counts at these positions in the tumor file, writes out a set of observables for use with bardCNV.

=cut

#------------
# Begin MAIN
#------------

process_commandline();
my $normal_file = $Opt{'normal_counts'};
my $tumor_file = $Opt{'tumor_counts'};

open NORM, $normal_file
    or die "Couldn\'t open $normal_file: $!\n";

my $rh_alt_allele = {};
my $rh_norm_count = {};
while (<NORM>) {
    if (/^(\S+)\s(\d+)\s(\S)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
        my ($chr, $pos, $ref, $A_count, $T_count, $G_count, $C_count, $a_count, $t_count, $g_count, $c_count) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
        next if ($ref !~ /[ATGC]/);
        my %base_counts = ();
        $base_counts{'A'} = $A_count + $a_count;
        $base_counts{'T'} = $T_count + $t_count;
        $base_counts{'C'} = $C_count + $c_count;
        $base_counts{'G'} = $G_count + $g_count;
        $rh_norm_count->{$chr}->{$pos} = $base_counts{'A'} + $base_counts{'T'} + $base_counts{'G'} + $base_counts{'C'};

        my $ref_count = $base_counts{$ref};
        delete $base_counts{$ref};

        my @sorted_bases = sort {$base_counts{$b} <=> $base_counts{$a}} keys %base_counts;
        my $alt_count = $base_counts{$sorted_bases[0]};
        $rh_alt_allele->{$chr}->{$pos} = $sorted_bases[0];
    }
}
close NORM;

open TUMOR, $tumor_file
    or die "Couldn\'t open $tumor_file: $!\n";

my $rh_alt_ratio = {};
my $rh_tumor_count = {};
while (<TUMOR>) {
    if (/^(\S+)\s(\d+)\s(\S)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
        my ($chr, $pos, $ref, $A_count, $T_count, $G_count, $C_count, $a_count, $t_count, $g_count, $c_count) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
        next if ($ref !~ /[ATGC]/);
        my %base_counts = ();
        $base_counts{'A'} = $A_count + $a_count;
        $base_counts{'T'} = $T_count + $t_count;
        $base_counts{'C'} = $C_count + $c_count;
        $base_counts{'G'} = $G_count + $g_count;
        $rh_tumor_count->{$chr}->{$pos} = $base_counts{'A'} + $base_counts{'T'} + $base_counts{'G'} + $base_counts{'C'};

        my $ref_count = $base_counts{$ref};
        my $alt_count;
        if ($rh_alt_allele->{$chr} && $rh_alt_allele->{$chr}->{$pos}) {
            my $alt_allele = $rh_alt_allele->{$chr}->{$pos};
            $alt_count = $base_counts{$alt_allele};
        }
        else {
            next;
        }

        $rh_alt_ratio->{$chr}->{$pos} = ($ref_count + $alt_count) ? $alt_count/($ref_count + $alt_count) : 0.5;
    }

}
close TUMOR;

foreach my $chr (keys %{$rh_norm_count}) {
    foreach my $pos (sort {$a <=> $b} keys %{$rh_norm_count->{$chr}}) {
        my $normal_reads = $rh_norm_count->{$chr}->{$pos};
        my $tumor_reads = ($rh_tumor_count->{$chr} && $rh_tumor_count->{$chr}->{$pos}) ? $rh_tumor_count->{$chr}->{$pos} : 0;
        next if ($normal_reads + $tumor_reads < $Opt{'mindepth'});
        my $rdr = $tumor_reads/$normal_reads;
        my $alt_freq = ($rh_alt_ratio->{$chr} && defined($rh_alt_ratio->{$chr}->{$pos})) ? $rh_alt_ratio->{$chr}->{$pos} : 0.5;
        print "$chr\t$pos\t$rdr\t$alt_freq\t$tumor_reads\t$normal_reads\n";
    }
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( mindepth => 0 );
    GetOptions(
        \%Opt, qw(
            manual help+ version normal_counts=s tumor_counts=s
            mindepth=i verbose 
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "make_observables.pl, ", q$Revision:$, "\n";
	}

}

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
