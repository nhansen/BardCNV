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

compare_cnv_beds.pl - script to read in two pseudo-BED files with CNV predictions and write out statistics on copy number agreement in regions covered in both.

=head1 SYNOPSIS

  compare_cnv_beds.pl <bedfile #1> <bedfile #2>

=head1 DESCRIPTION

This script reads in two files with tab-delimited columns (chr, start, end, minor copy number, total copy number, quality score), and reports statistics regarding the agreement of copy number predictions in regions covered by both files.

=cut

#------------
# Begin MAIN
#------------

process_commandline();

my $bedfile1 = $ARGV[0];
my $bedfile2 = $ARGV[1];

my $rh_regions = {}; # ref to hash on chromosome, values are refs to lists of hashes with start, end, minor, total, score

open BED1, $bedfile1
    or die "Couldn\'t open $bedfile1: $!\n";

while (<BED1>) {
    if (/^(\S+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)\s(\d+)$/) {
         my ($chr, $start, $end, $minor, $total, $score) = ($1, $2, $3, $4, $5, $6);
         if (!$Opt{'minscore'} || $score >= $Opt{'minscore'}) {
             push @{$rh_regions->{$chr}}, {'start' => $start,
                                           'end' => $end,
                                           'minor' => $minor,
                                           'total' => $total,
                                           'score' => $score };
         }
    }
}

close BED1;

my $rh_overlap = {};
my $total_overlap = 0;

open BED2, $bedfile2
    or die "Couldn\'t open $bedfile2: $!\n";

while (<BED2>) {
    if (/^(\S+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)\s(\d+)$/) {
         my ($chr, $start, $end, $minor, $total, $score) = ($1, $2, $3, $4, $5, $6);

         next if (!$rh_regions->{$chr});
         my @regions_to_remove = ();
         if (!$Opt{'minscore'} || $score >= $Opt{'minscore'}) {
             # search chromosome's regions for overlap:

             my $remaining_regions = $#{$rh_regions->{$chr}};
             for (my $reg_index = 0; $reg_index <= $#{$rh_regions->{$chr}}; $reg_index++) {
                 my $rh_reg_info = $rh_regions->{$chr}->[$reg_index];
                 my $reg_start = $rh_reg_info->{'start'}; 
                 my $reg_end = $rh_reg_info->{'end'}; 
                 my $reg_minor = $rh_reg_info->{'minor'}; 
                 my $reg_total = $rh_reg_info->{'total'}; 
                 my $reg_score = $rh_reg_info->{'score'}; 

                 if ($start > $reg_end) { # can remove region
                     push @regions_to_remove, $reg_index;
                 }
                 elsif ($end > $reg_start) { # this bed1 region overlaps with our bed2 region
                     if (!$Opt{'minscore'} || $score >= $Opt{'minscore'}) {
                         my $overlap_start = ($reg_start >= $start) ? $reg_start : $start;
                         my $overlap_end = ($reg_end >= $end) ? $end : $reg_end;
                         $total_overlap += $overlap_end - $overlap_start;
                         my $total_diff = $total - $reg_total;
                         $rh_overlap->{$total_diff} += $overlap_end - $overlap_start;

                         if (defined($Opt{'report'}) && abs($total_diff) >= $Opt{'report'}) {
                             my $overlap_length = $overlap_end - $overlap_start;
                             print "DISC\t$chr\t$overlap_start\t$overlap_end\t$overlap_length\t$reg_minor\t$reg_total\t$reg_score\t$minor\t$total\t$score\n";
                         }
                     }
                     if ($reg_end > $end) { # no more regions to search
                         last;
                     }
                 }
             }
             my $splice_offset = 0;
             foreach my $index (@regions_to_remove) {
                 splice(@{$rh_regions->{$chr}}, $index - $splice_offset, 1);
                 $splice_offset++;
             }
         }
    }
}
close BED2;

foreach my $diff (sort {$a <=> $b} keys %{$rh_overlap}) {
    my $overlap_bases = $rh_overlap->{$diff};
    my $overlap_percent = $overlap_bases/$total_overlap;
    print "$diff\t$overlap_bases\t$total_overlap\t$overlap_percent\n";
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( minscore => 0 );
    GetOptions(
        \%Opt, qw(
            manual help+ version minscore=i report=i verbose 
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "compare_cnv_beds.pl, ", q$Revision:$, "\n";
	}

}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--minscore>

Specify the minimum score a region needs to be counted as "covered".

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
