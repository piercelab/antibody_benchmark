#!/usr/bin/perl

use strict;

my $file = $ARGV[0];
if ($file eq "") { die("usage: get_cases.pl antibody_antigen_triplet_list\n"); } 
#"bound_unbound_pdbs.092020.txt";

my $pdb_res_file = "/Users/bpierce/databases/pdb_resolutions_8_2020.txt";
my $current_cases = "bm5.5.txt";

my @curr_bm = ();

open(FL, $file) || die("unable to open cases file: $file\n");
my @file_lines = <FL>;
close(FL);
chomp(@file_lines);

my %possible_cases = ();
foreach my $line (@file_lines)
{
    my @fields = split("\t", $line);
    if (($fields[1] ne "") && ($fields[2] ne "")) { $possible_cases{$fields[0]} = "$fields[1]\t$fields[2]"; }
}

open(RES, $pdb_res_file) || die("unable to open file: $pdb_res_file\n");
my @res_lines = <RES>;
close(RES);

# parse those PDB resolutions
my %pdbres = ();
my $no_res_yet = 1;
for (my $i = 0; $i < @res_lines; $i++)
{
    if (substr($res_lines[$i], 0, 6) eq "IDCODE") { $i += 2; $no_res_yet = 0; }
    if ($no_res_yet == 1) { next; }
    (my $pdb, my $stuff, my $res) = split(" ", $res_lines[$i]);
    $pdbres{lc($pdb)} = $res;
}

open(BM5, $current_cases) || die("unable to open file: $current_cases\n");
my @bm5_lines = <BM5>;
close(BM5);
chomp(@bm5_lines);

my $bm5_matches = 0;
foreach my $line (@bm5_lines)
{
    if ($possible_cases{$line} ne "") { $bm5_matches++; push @curr_bm, $line; }
}
print $bm5_matches . " out of " . scalar(@bm5_lines) . " BM5 matches\n";

# trim those cases down
my @all_comps = keys %possible_cases;
foreach my $comp (@all_comps)
{
    (my $mabs, my $ags) = split("\t", $possible_cases{$comp});
    my @mab_comps = split(" ", $mabs);
    my $top_mab = $mab_comps[0];
    my $top_mab_res = $pdbres{$mab_comps[0]};
    for (my $i = 1; $i < @mab_comps; $i++) { if ($pdbres{$mab_comps[$i]} < $top_mab_res) { $top_mab = $mab_comps[$i]; $top_mab_res = $pdbres{$mab_comps[$i]}; } }
    my @ag_comps = split(" ", $ags);
    my $top_ag = $ag_comps[0];
    my $top_ag_res = $pdbres{$ag_comps[0]};
    for (my $i = 1; $i < @ag_comps; $i++) { if ($pdbres{$ag_comps[$i]} < $top_ag_res) { $top_ag = $ag_comps[$i]; $top_ag_res = $pdbres{$ag_comps[$i]}; } }
    $possible_cases{$comp} = $top_mab . "\t" . $top_ag;
}


my $new_cases = 0;
my @new_case_comps = ();
foreach my $comp (@all_comps)
{
    my $match_found = 0;
    for (my $i = 0; $i < @curr_bm; $i++)
    {
	#if (($comp eq $curr_bm[$i]) || ($possible_cases{$comp} eq $possible_cases{$curr_bm[$i]})) { $match_found = 1; last; }
	if ($comp eq $curr_bm[$i]) { $match_found = 1; last; }
    }
    if ($match_found == 0) 
    { 
	push @curr_bm, $comp; 
	push @new_case_comps, $comp; 
	$new_cases++; 
    }
}
print $new_cases . " new cases found!\n";
foreach my $case (@new_case_comps) { print $case . "\t" . $possible_cases{$case} . "\n"; }
