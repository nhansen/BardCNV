#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw/ tempdir tempfile /;
use File::Copy;
use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr( q$Rev: 0$, 4 ) / 10000;

=head1 NAME

bardCNV.pl - given a file of heterozygous positions (in a normal diploid sample), set the parameters and solve for the states of a hidden Markov model to predict copy number state in a matched tumor sample.

=head1 SYNOPSIS

  bardCNV.pl --normalbam <bam file from normal sample> --tumorbam <bam file from tumor> --ref <reference fasta> --hetfile <BED formatted file of positions heterozygous in the normal>
  bardCNV.pl --trio --mombam <bam file from mother's sample> --dadbam <bam file from father's sample> --childbam <bam file from child>

=head1 DESCRIPTION

Using B-allele frequencies within sequence from a tumor sample, and the depth of coverage of the tumor compared to normal sequencing read depth, bardCNV.pl predicts copy number (including the copy counts of both alleles) using an HMM trained on these two features.  The program requires a list of high-confidence sites which are believed to be diploid and heterozygous in the normal sample.  For these sites, the B-allele frequency in the tumor sample helps bardCNV to pin down the ratio of the major to minor copy number allele, and thus helps to determine the overall ploidy, and optionally, the purity level of the tumor sample.

=cut

#------------
# Begin MAIN
#------------

$ENV{'PATH'} = "/home/nhansen/projects/bard_binomial/bard/scripts:/home/nhansen/projects/bard_binomial/bard/c:/home/nhansen/projects/bamcounts:$ENV{PATH}";

my $plotstates_rlib = "/home/nhansen/projects/bard_binomial/bard/R/plot_states.R";

process_commandline();

my $workdir = create_working_directory();
my ($train_observable_file, $state_observable_file) = create_observable_files($workdir);
my $best_param_file = train_model($workdir, $train_observable_file);
my $state_file = run_viterbi($workdir, $state_observable_file, $best_param_file);
my $vcf_file = write_vcf($workdir, $state_observable_file, $state_file);
my $stats_file = calc_stats_file($workdir, $state_file);

if ($Opt{'plots'}) {
    my $plotdir = draw_plots($workdir, $best_param_file, $state_file);
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( mindepth => 0 , maxcopies => 8, transprob => 0.0001, countopts => '',
             sigratio => 3.00, contam => 0.01, sigpi => 0.15, maxviterbicopies => 24,
             mincontam => 0.01, maxcontam => 0.5 );
    GetOptions(
        \%Opt, qw(
            manual help+ version normalbam=s tumorbam=s ref=s
            outdir=s hetfile=s mindepth=i maxcopies=i 
            transprob=f countopts=s trio fixedtrans
            mombam=s dadbam=s childbam=s male diploid
            sigratio=f contam=f optcontam mincontam=i
            maxcontam=i sigpi=f skipcounts skipobs 
            skiptrain skipstates skipvcf verbose
            omittrainentries=s plots sge vcf!
		  )
	) || pod2usage(0);
    if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
    if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
    if ( $Opt{version} ) {
        die "bardCNV.pl, ", q$Revision:$, "\n";
    }
    if ($Opt{'trio'}) {
        if (!$Opt{'mombam'} || !$Opt{'dadbam'} || !$Opt{'childbam'}) {
            die "Options --mombam, --dadbam, --childbam, and --ref are required with --trio option!\n";
        }
    }
    elsif (!$Opt{'normalbam'} || !$Opt{'tumorbam'} || !$Opt{'ref'} || !$Opt{'hetfile'}) {
        die "Options --normalbam, --tumorbam, --hetfile, and --ref are required!\n";
    }
    if ($Opt{'optcontam'}) {
        if ($Opt{'mincontam'} < 0.01) {
           $Opt{'mincontam'} = 0.01;
        }
        if ($Opt{'maxcontam'} > 0.99) {
           $Opt{'maxcontam'} = 0.99;
        }
    }
}

sub create_working_directory {

    my $workdir;
    if ($Opt{outdir}) {
        if (!(-e $Opt{outdir})) {
            print STDERR "Attempting to create $Opt{'outdir'}\n" if ($Opt{'verbose'});
            mkdir $Opt{'outdir'};
        }
        $workdir = $Opt{'outdir'};
    }
    else {
        $workdir = tempdir("run_bardcnv_XXXXXX", DIR => '.');
    }

    if ($workdir !~ /^\//) {
        my $pwd = `pwd`;
        chomp $pwd;
        $workdir = "$pwd/$workdir";
    }

    return $workdir;
}

sub create_observable_files {
    my $workdir = shift;

    my $train_obs_file = "$workdir/train_obs.txt";
    my $viterbi_obs_file = "$workdir/viterbi_obs.txt";

    if ($Opt{'skipobs'}) {
        if (!(-e $train_obs_file) || !(-e $viterbi_obs_file)) {
            die "Files $train_obs_file and $viterbi_obs_file must exist when option --skipobs is specified!\n";
        }
        return ($train_obs_file, $viterbi_obs_file);
    }

    if ($Opt{'skipcounts'}) {
        if (!(-e "$workdir/normalbam.counts") || !(-e "$workdir/tumorbam.counts")) {
            die "Files $workdir/normalbam.counts and $workdir/tumorbam.counts must exist when option --skipcounts is specified!\n";
        }
    }
    else { # generate counts files:
        my @count_commands = ();
        my @bam_files = ($Opt{'trio'}) ? ('mombam', 'dadbam', 'childbam') : ('normalbam', 'tumorbam');
        foreach my $fileopt (@bam_files) {
            my $cmd = "bamcounts -fasta $Opt{ref} -bam $Opt{$fileopt} -bedfile $Opt{hetfile} $Opt{countopts} > $workdir/$fileopt.counts";
            push @count_commands, $cmd;
        }
        run_commands(\@count_commands, '-o' => "$workdir",
                           '-e' => "$workdir"); 
    }
    
    my $rh_norm_count = {};
    my $rh_tumor_count = {};
    my $rh_alt_count = {}; # the number of reads displaying the alternate allele in the tumor
    my $rh_alt_allele = {}; # this will be the non-reference allele which is present in the most reads in the normal

    my @normal_count_files = ($Opt{'trio'}) ? ("$workdir/mombam.counts", "$workdir/dadbam.counts") :
                                              ("$workdir/normalbam.counts");
    foreach my $normal_file (@normal_count_files) {
        open NORM, "$normal_file"
            or die "Couldn\'t open $normal_file: $!\n";
        
        while (<NORM>) {
            if (/^(\S+)\s(\d+)\s(\S)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
                my ($chr, $pos, $ref, $A_count, $T_count, $G_count, $C_count, $a_count, $t_count, $g_count, $c_count) = ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11);
                next if ($ref !~ /[ATGC]/);
                my %base_counts = ();
                $base_counts{'A'} = $A_count + $a_count;
                $base_counts{'T'} = $T_count + $t_count;
                $base_counts{'C'} = $C_count + $c_count;
                $base_counts{'G'} = $G_count + $g_count;

                # add in coverage, multiplying by two if this is a dad on the X chrom in a trio:
                my $multiplier = 1;
                if (($Opt{'trio'}) && ($normal_file eq "$workdir/dadbam.counts") && 
                                 ($chr =~ /x/i) && (nonpar_xpos($pos))) {
                    $multiplier = 2;
                }

                $rh_norm_count->{$chr}->{$pos} += $multiplier * ($base_counts{'A'} + $base_counts{'T'} + 
                                                        $base_counts{'G'} + $base_counts{'C'});
        
                if (!$Opt{'trio'}) {
                    delete $base_counts{$ref}; # for trios, only one sample will display alternate allele, so don't want to record wrong alternate allele
                }
        
                my @sorted_bases = sort {$base_counts{$b} <=> $base_counts{$a}} keys %base_counts;
                my $alt_count = $base_counts{$sorted_bases[0]};
                if ($sorted_bases[0] ne $ref) {
                    $rh_alt_allele->{$chr}->{$pos} = $sorted_bases[0];
                }
            }
        }
        close NORM;
    }
   
    my $tumorfile = ($Opt{'trio'}) ? "$workdir/childbam.counts" : "$workdir/tumorbam.counts"; 
    open TUMOR, $tumorfile
        or die "Couldn\'t open $tumorfile: $!\n";
    
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
    
            $rh_alt_count->{$chr}->{$pos} = $alt_count;
        }
    }
    close TUMOR;

    my @omit_entries = ($Opt{'omittrainentries'}) ? split ',', $Opt{'omittrainentries'} : ();
  
    open STATEOBS, ">$viterbi_obs_file"
        or die "Couldn\'t open $viterbi_obs_file for writing: $!\n";

    open TRAINOBS, ">$train_obs_file"
        or die "Couldn\'t open $train_obs_file for writing: $!\n";

    foreach my $chr (keys %{$rh_norm_count}) {
        foreach my $pos (sort {$a <=> $b} keys %{$rh_norm_count->{$chr}}) {
            my $normal_reads = $rh_norm_count->{$chr}->{$pos};
            my $tumor_reads = ($rh_tumor_count->{$chr} && $rh_tumor_count->{$chr}->{$pos}) ? $rh_tumor_count->{$chr}->{$pos} : 0;
            next if ($normal_reads + $tumor_reads < $Opt{'mindepth'});
            next if (!$normal_reads);
            my $alt_count = ($rh_alt_count->{$chr} && defined($rh_alt_count->{$chr}->{$pos})) ? $rh_alt_count->{$chr}->{$pos} : 0;
            if ((!@omit_entries || !grep {$_ eq $chr} @omit_entries) && ($chr !~ /X/ || !$Opt{'male'} )) {
                print TRAINOBS "$chr\t$pos\t$normal_reads\t$tumor_reads\t$alt_count\n";
            }
            print STATEOBS "$chr\t$pos\t$normal_reads\t$tumor_reads\t$alt_count\n";
        }
    }
    close STATEOBS;
    close TRAINOBS;

    return ($train_obs_file, $viterbi_obs_file);
}

sub train_model {
    my $workdir = shift;
    my $obs_file = shift;

    # calculate mean of tumor read count distribution:
    open OBS, $obs_file
        or die "Couldn\'t open $obs_file: $!\n";
    my $total_ratio = 0.0;
    my $datapoints = 0;
    my $total_normal_depth = 0;
    my $total_tumor_depth = 0;
    while (<OBS>) {
        if (/^(\S+)\s(\d+)\s(\d+)\s(\d+)/) { # read depth ratio is fourth divided by third column
            my ($normal_depth, $tumor_depth) = ($3, $4);
            if (!$normal_depth) {
                die "Normal depth of coverage is $normal_depth at $1:$2--this shouldn\'t happen!\n";
            }
           
            my $ratio = $tumor_depth/$normal_depth; 
            $total_ratio += $ratio;
            $datapoints++;
            $total_normal_depth += $normal_depth;
            $total_tumor_depth += $tumor_depth;
        }
    }
    close OBS;
    my $ratio_avg = $total_ratio/$datapoints;
    print STDERR "Mean read depth ratio: $ratio_avg\n";
    my @muratio_vals = ($ratio_avg/4.0, $ratio_avg/3.0, $ratio_avg/2.0, 2.0*$ratio_avg/3.0, $ratio_avg, 2.0*$ratio_avg);
    my @contam_vals = ();
    if ($Opt{'optcontam'}) {
        for (my $thiscontam = $Opt{'mincontam'}; $thiscontam <= $Opt{'maxcontam'}; $thiscontam+= 0.01) {
            push @contam_vals, $thiscontam;
        }
    }
    else {
        push @contam_vals, $Opt{'contam'};
    }

    my @train_commands = ();
    my @output_files = ();
    foreach my $muratio (@muratio_vals) {
        my $muratio_string = sprintf("%5.3f", $muratio);
        foreach my $contamval (@contam_vals) {
            if (!$Opt{'skiptrain'}) {
                if (!(-e "$workdir/train_$muratio_string\_$contamval")) {
                    mkdir "$workdir/train_$muratio_string\_$contamval"
                        or die "Couldn\'t create $workdir/train_$muratio_string\_$contamval: $!\n";
                }
                write_train_start_file("$workdir/train_$muratio_string\_$contamval/train_start.txt", $muratio, $contamval); 
                my $maxratio_opt = $ratio_avg * 4.0;
                my $fix_opt = ($Opt{'fixedtrans'}) ? '-fixtrans' : '';
                my $verbose_opt = ($Opt{'verbose'}) ? '-verbose' : '';
                my $train_cmd = "bardcnv baumwelch -derivatives -obsfile $obs_file -modelfile $workdir/train_$muratio_string\_$contamval/train_start.txt $fix_opt $verbose_opt -maxratio $maxratio_opt > $workdir/train_$muratio_string\_$contamval/train.out 2>$workdir/train_$muratio_string\_$contamval/train.err";
                push @train_commands, $train_cmd;
            }
            push @output_files, "$workdir/train_$muratio_string\_$contamval/train.out";
        }
    }
    if (!$Opt{'skiptrain'}) {
        run_commands(\@train_commands, -mem_free => '80G', '-o' => $workdir, '-e' => $workdir);
    }

    my $best_model_file;
    my $bestlogp;
    my $diploid_model_file;
    my $most_diploid_ploidy;
    foreach my $outputfile (@output_files) {
        my $error_file = $outputfile;
        $error_file =~ s/\.out/.err/;
        my $finished = `grep 'Finished' $error_file`;
        if (!$finished) {
            print STDERR "Skipping training file $error_file--did not complete successfully\n";
            next;
        }
        my $logp = `grep 'MLL difference' $error_file | tail -1 | awk '{print \$5}' | sed 's/)//'`;
        chomp $logp;
        my $mu = `grep 'mu=' $outputfile | awk '{print \$2}'`;
        chomp $mu;
        my $contam = `grep 'contam=' $outputfile | awk '{print \$2}'`;
        chomp $contam;
        print STDERR "$error_file\t$mu\t$contam\t$logp\n";

        if (($logp =~ /[-\d\.]+/) && (!(defined($bestlogp)) || ($logp > $bestlogp))) {
            $best_model_file = $outputfile;
            $bestlogp = $logp;
        }

        if (($Opt{'diploid'}) && $mu =~ /^[.\d]+$/ && $contam =~ /^[.\d]+$/) {
            my $ploidy = (2.0*$total_tumor_depth/$total_normal_depth/$mu - $contam*2.0)/(1.0 - $contam);
            if (!(defined($most_diploid_ploidy)) || abs($ploidy - 2.0) < abs($most_diploid_ploidy - 2.0)) {
                $most_diploid_ploidy = $ploidy;
                $diploid_model_file = $outputfile;
            }
        }
    }
    if ($best_model_file) {
        copy($best_model_file, "$workdir/best_model_file.txt");
    }
    else {
        die "Unable to train model--no best mu_ratio value.\n";
    }

    if ($diploid_model_file) {
        copy($diploid_model_file, "$workdir/diploid_model_file.txt");
    }
    else {
        die "Unable to create a diploid model file--no training produced a close-to-diploid model.\n";
    }

    return ($Opt{'diploid'}) ? "$workdir/diploid_model_file.txt" : "$workdir/best_model_file.txt";
}

sub write_train_start_file {
    my $startfile = shift;
    my $muratio = shift;
    my $contamval = shift;

    my $maxcopies = $Opt{'maxcopies'};
    my $total_states = 0;
    for (my $copies = 0; $copies <= $maxcopies; $copies++) {
        $total_states += int($copies/2.0 + 1);
    }

    # allocation of initial state probabilities--probably not very important when considering an entire exome or genome...
    my $largeprob = 0.8; # probability of normal, diploid state
    my $smallprob = (1.0 - $largeprob)/($total_states - 1.0); # probability of any other state

    open START, ">$startfile"
        or die "Couldn\'t open $startfile for writing: $!\n";

    print START "N= $total_states\n"; 
    for (my $copies = 0; $copies <= $maxcopies; $copies++) {
        for (my $minor = 0; $minor < int($copies/2.0 + 1); $minor++) {
            if ($minor == 1 && $copies == 2) {
                print START "$minor $copies $largeprob\n";
            }
            else {
                print START "$minor $copies $smallprob\n";
            }
        }
    }

    print START "Transprob= $Opt{'transprob'}\n";
    print START "mu= $muratio\n";
    print START "sigma= $Opt{'sigratio'}\n";
    print START "rho_contam= $contamval\n";

    close START;
}

sub run_viterbi {
    my $workdir = shift;
    my $obs_file = shift;
    my $param_file = shift;

    if ($Opt{'skipstates'}) { # verify state file exists and return its path:
        if (-e ("$workdir/states.out")) {
            return "$workdir/states.out";
        }
        else {
            die "Option --skipstates invalid--file $workdir/states.out doesn\'t exist!\n";
        }
    }

    # add high copy number states to model:
    my $orig_states;
    open TRAIN, "$param_file"
        or die "Couldn\'t open $param_file to add high copy number states: $!\n";

    my $viterbi_param_file = $param_file;
    $viterbi_param_file =~ s/\.txt$//;
    $viterbi_param_file .= '.highcopy.txt';

    open NEWTRAIN, ">$viterbi_param_file"
        or die "Couldn\'t open $viterbi_param_file for writing: $!\n";

    print STDERR "Opening file $viterbi_param_file for writing high copy number model params.\n";

    my @state_lines = ();
    while (<TRAIN>) {
        if (/^N= (\d+)$/) {
            $orig_states = $1;
            for (my $state = 0; $state < $orig_states; $state++) {
                my $next_line = <TRAIN>;
                push @state_lines, $next_line;
            }
            my $last_line = $state_lines[$#state_lines];
            if ($state_lines[$#state_lines] =~ /^(\d+)\s(\d+)/) {
                my $train_max_states =  $2;
                for (my $copies = $train_max_states + 1; $copies <= $Opt{'maxviterbicopies'}; $copies++) {
                    for (my $minor = 0; $minor < int($copies/2.0 + 1); $minor++) {
                        push @state_lines, "$minor $copies 0.0\n";
                    }
                }
                my $new_states = @state_lines;
                print NEWTRAIN "N= $new_states\n";
                foreach my $line (@state_lines) {
                    print NEWTRAIN $line;
                }
            }
            else {
                close TRAIN;
                die "Misformatted training file $param_file\n";
            }
        }
        else {
            print NEWTRAIN $_;
        }
    }
    close NEWTRAIN;
  
    my $cmd = "bardcnv viterbi -obsfile $obs_file -modelfile $viterbi_param_file > $workdir/states.out 2>$workdir/states.err";
    run_command($cmd, 'mem_free' => '40G', '-o' => $workdir, '-e' => $workdir);

    return "$workdir/states.out";
}

sub write_vcf {
    my $workdir = shift;
    my $obs_file = shift;
    my $state_file = shift;

    my $vcf_file = "$workdir/states.vcf.gz";
    if ($Opt{'skipvcf'}) {
        return $vcf_file;
    }

    my $cmd = "bard2vcf.pl --obs $obs_file --states $state_file --ref $Opt{'ref'} --sample1 NORMAL --sample2 TUMOR --outfile $vcf_file";
    run_command($cmd, '-o' => $workdir, '-e' => $workdir);

    return $vcf_file;
}

sub calc_stats_file {
    my $workdir = shift;
    my $state_file = shift;

    my $stats_file = "$workdir/qual_stats.txt";

    my $rh_qual_stats = {};
    my %all_chroms = ();
    open STATES, $state_file
        or die "Couldn\'t open $state_file: $!\n";
    while (<STATES>) {
        if (/^(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my $chr = $1;
            my $score = $8;
            $all_chroms{$chr} = 1;
            $rh_qual_stats->{$chr}->{'total'} = 0 if (! (defined ($rh_qual_stats->{$chr}->{'total'})));
            $rh_qual_stats->{$chr}->{'60'} = 0 if (! (defined ($rh_qual_stats->{$chr}->{'60'})));
            $rh_qual_stats->{$chr}->{'90'} = 0 if (! (defined ($rh_qual_stats->{$chr}->{'90'})));
            $rh_qual_stats->{$chr}->{'total'}++;
            $rh_qual_stats->{$chr}->{'60'}++ if ($score>=60);;
            $rh_qual_stats->{$chr}->{'90'}++ if ($score>=90);;
        }
    }
    close STATES;

    open STATS, ">$stats_file"
        or die "Couldn\'t open $stats_file for writing: $!\n";

    print STATS "Chrom\tTotalSites\tPercent60orAbove\tPercent90orAbove\n";
    foreach my $chr (sort {$rh_qual_stats->{$b}->{'60'}/$rh_qual_stats->{$b}->{'total'} <=> 
                    $rh_qual_stats->{$a}->{'60'}/$rh_qual_stats->{$a}->{'total'}} keys %all_chroms) {
        my $total_points = $rh_qual_stats->{$chr}->{'total'};
        my $perc60 = int(1000*$rh_qual_stats->{$chr}->{'60'}/$total_points)/10;
        my $perc90 = int(1000*$rh_qual_stats->{$chr}->{'90'}/$total_points)/10;

        print STATS "$chr\t$total_points\t$perc60\t$perc90\n";
    }
    close STATS;

    return $stats_file;
}

sub draw_plots {
    my $workdir = shift;
    my $param_file = shift;
    my $state_file = shift;

    mkdir("$workdir/plots");

    my $muratio = `grep 'mu=' $param_file | awk '{print \$2}'`;
    chomp $muratio;
    my $contam = `grep 'rho_contam=' $param_file | awk '{print \$2}'`;
    chomp $contam;

    my $r_cmd_file = "$workdir/plots/plotstategraphs.R";
    open COM, ">$r_cmd_file"
        or die "Couldn\'t open $r_cmd_file for writing: $!\n";

    print COM <<"DOC";
source("$plotstates_rlib");

states <- read_states("$state_file");

states60 <- states[states\$score >= 60, ];
chroms <- unique(states\$chr);

for (thischr in chroms) {
    file = paste("$workdir/plots/", thischr, ".png", sep="");
    png(file);
    plot_states(states, thischr, muratio=$muratio, contam=$contam);
    dev.off();
}

chromshq <- unique(states60\$chr);

for (thischr in chromshq) {
    file = paste("$workdir/plots/", thischr, ".hq.png", sep="");
    png(file);
    plot_states(states60, thischr, muratio=$muratio, contam=$contam);
    dev.off();
}

DOC

    my $cmd = "R --file=$r_cmd_file";
    run_command($cmd, '-o' => $workdir, '-e' => $workdir);
}

sub run_command {
    my $cmd = shift;
    my %params = @_;

    my $rh_jobinfo;
    if ($Opt{'sge'}) {
        my $mem_required = ($params{'-mem_free'}) ? $params{'-mem_free'} : '1G';
        delete $params{'-mem_free'} if ($params{'-mem_free'});
        $rh_jobinfo = sge_submit($cmd, 'mem_required' => $mem_required, %params);
    }
    else {
        system($cmd);
    }

    if ($Opt{'sge'}) { # wait for submitted jobs to finish
        system("qrsh -now no -hold_jid $rh_jobinfo->{'jobid'} -b y /bin/true");
        unlink $rh_jobinfo->{cmdfile};
    }
}

sub run_commands {
    my $ra_cmds = shift;
    my %params = @_;

    my $mem_required = ($params{'-mem_free'}) ? $params{'-mem_free'} : '1G';
    delete $params{'-mem_free'} if ($params{'-mem_free'});
 
    my @sge_job_ids = ();
    my @cmd_files = ();
    foreach my $cmd (@{$ra_cmds}) { 
        if ($Opt{'sge'}) {
            my $rh_jobinfo = sge_submit($cmd, 'mem_required' => $mem_required, %params);
    
            push @cmd_files, $rh_jobinfo->{cmdfile};
            if ($rh_jobinfo->{jobid} =~ /^(\d+)$/) {
                push @sge_job_ids, $rh_jobinfo->{jobid};
            }
            else {
                die "Something went wrong submitting job $cmd\n";
            }
        }
        else {
            system($cmd);
        }
    }
    if ($Opt{'sge'}) { # wait for submitted jobs to finish
        my $hold_string = join ',', @sge_job_ids;
        system("qrsh -now no -hold_jid $hold_string -b y /bin/true");
        foreach my $cmd_file (@cmd_files) {
            unlink $cmd_file;
        }
    }
}

sub sge_submit {
    my $cmd = shift;
    my %params = @_;

    my $mem_required = $params{'mem_required'} || '1G';

    my $jobname = $cmd;
    $jobname =~ s/\s.*$//;
    my $memlimit = $mem_required;
    $memlimit =~ s/G/000000/;
    my ($fh, $tempfile) = tempfile("tmpXXXXXX", DIR => '.');
    print $fh "#!/bin/bash\nulimit $memlimit\n$cmd\n";
    $fh->close();
    my $qsub_cmd = "qsub -terse -l mem_free=$mem_required -N $jobname";
    if ($params{'-o'}) {
        $qsub_cmd .= " -o $params{'-o'}";
    }
    if ($params{'-e'}) {
        $qsub_cmd .= " -e $params{'-e'}";
    }
    $qsub_cmd .= " $tempfile | ";
    open QSUB, "$qsub_cmd"
        or die "Couldn\'t open qsub for $cmd\n";
    my $job_id = <QSUB>;
    close QSUB;
    chomp $job_id;

    return {'cmdfile' => $tempfile, 'jobid' => $job_id};
}

sub nonpar_xpos {
    my $pos = shift;

    if ($Opt{'ref'} =~ /hg18/i) {
        if ($pos >= 2709521 && $pos <= 154584238) {
            return 1;
        }
    }
    elsif ($Opt{'ref'} =~ /hg19/i) {
        if ($pos >= 2699521 && $pos <= 154931044) {
            return 1;
        }
    }
    else {
        die "Unable to determine non-PAR regions of X-chromosome: Can\'t determine build from reference fasta name!\n";
    }

    return 0;
}


__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--normalbam>

Specify the path of a BAM file containing the aligned sequence data from the normal sample (required).

=item B<--tumorbam>

Specify the path of a BAM file containing the aligned sequence data from the tumor sample (required).

=item B<--ref>

Specify the path of the reference fasta file to which reads from the BAM files are aligned (required).

=item B<--mindepth>

Specify the minimum total depth of coverage required to included a heterozygous site in the observables written to the observable file (default 0).

=item B<--maxcopies>

Specify the maximum total number of copies to use in the training step (limits memory and CPU usage) (default 12, maximum 16).

=item B<--transprob>

Specify the starting value for the transition probability between different states (which is constrained to be the same for all inter-state transitions (default 0.0001).  To fix the transition probability so that it is not optimized using Baum-Welch, use the "--fixedtrans" option.

=item B<--fixedtrans>

Do not optimize the state transition probabilities using Baum-Welch.  Often, allowing the transition probability to vary will result in overfitting of the model with an overly-large transition probability, allowing the introduction of small artifactual copy number changes to increase the likelihood of data that doesn't conform to the overly-simplistic Gaussian model.

=item B<--sigratio>

Specify the starting value for the standard deviation of the (Gaussian distributed) ratio of depth of coverage in the tumor sequence to depth of coverage in the normal sequence (default 3.00).

=item B<--sigpi>

Specify the starting value for the standard deviation of the (Gaussian distributed) proportion of alternate allele in the tumor sequence (default 0.15).

=item B<--contam>

Specify the percentage contamination of the tumor tissue which is estimated to actually be unaffected (default 0.01).

=item B<--optcontam>

Train the model for a range of contamination values, then output states for the optimal model at the optimal contamination.

=item B<--mincontam>

If using --optcontam option, this is the minimum value of the contamination parameter to test (default 0.01).  Note: if any value less then 0.01 is specified, 0.01 will be used.

=item B<--maxcontam>

If using --optcontam option, this is the maximum value of the contamination parameter to test (default 0.50).  Note: if any value greater than 0.99 is specified, 0.99 will be used.

=item B<--male>

Omit any chromosome containing the character "X" in the training step (this is because sites, at least in the non-PAR regions of X for a male, cannot be heterozygous.

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
