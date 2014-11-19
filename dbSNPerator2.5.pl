#!/usr/bin/perl -w

#Cohn lab script for adding dbSNP data to seattleseq outputs.  Requires splitting the reference VCF with splitterlibrary.
#Version info: V2.5 - Changed reference data architecture to use a split library of chromosomes
#instead of whole chromosomes.  Increased speed significantly.  Added loading of whole library file into
#memory to further increase speed to Nimbus 2000 levels.
#Copyright 2014 Michael Weinstein, Daniel H. Cohn laboratory, UCLA.

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);
$|=1;


sub main(){

    my %opts;  #setup the hash for commandline options
    getopts('f:', \%opts);  #write the commandline options to the array
    unless(checkopts(\%opts)){  #check the options for validity
        usage();  #if they are not valid (function returns false), print the usage instructions
    }
    
    unless(checkintegrity()){ #runs the integrity check subroutine to make sure that the subvcf directory contains the appropriate files
        die "Subvcf directory integrity failed. Likely due to missing or modified files\.\n"; #quits if it does not
    }
    my $seqin = $opts{f};  #sets up the variable for the file to be annotated using the commandline option
    my $seqout = $seqin."\.SNPed\.txt";  #sets up the variable for the output file for annotated data
    open(SSI, $seqin) or die "Couldn't open file $seqin. \n"; #opens the input file or quits
    if (-e $seqout) {  #checks to see if the annotation output file already exists
        die "\nOutput file $seqout already appears to exist\.\n"; #quits if it does
    }
    if (-e $seqout."\.log") { #checks to see if the annotation log file already exists
        die "\nLog file $seqout\.log already appears to exist\.\n"; #quits if it does
    }
    
 

    open(OUTPUT, '>'.$seqout) or die "\nError opening output file $seqout\.\n";  #Creates the file for annotation output
    open(LOGFILE, '>'.$seqout."\.log") or die "\nError opening log file $seqout\.log\.\n"; #Creates the log file
    my $header = <SSI>;  #reads the first line of the file (which should contain header data)  #reads the first seattleseq output line as the header
    $header =~ s/\R//g;  #removes any linebreaks from the header (including the break at the end)
    print OUTPUT $header."\tdbSNP ID\tdbSNP 5\% in all\tdbSNP 5\% in any\tdbSNP Mutation\n";  #Writes the new header line to the annotation output file
    unless ($header =~ /chromosome\tposition/i){  #checks to see if the header looks anything like a seattleseq header
        die "File header for $seqin does not indicate chromosome and position data." #quits if it doesn't
    }
    my @header = split(/\t/, $header);  #splits the header on each tab into a separate variable in an array
    my $progress = 0; #starts the progress counter at zero
    my $lastmillimorgan = -1;  #initializes the millimorgan variable to -1 (an impossible value once it begins reading)
    my $snplibrary;  #initializes a variable to store the entirety of the current dbSNP vcf portion
    
    LINE: while (my $line = <SSI>){  #loop that goes through file while there is file to go through line by line
        
        #open(OUTPUT, '>>'.$seqout) or die "\nError opening output file $seqout\.\n";
        if ($line =~ /^#/) {  #checks to see if the line starts with a #, indicating a data or comment line
            print OUTPUT "$line"; #if so, prints it directly to the output file unmodified
            next LINE; #and moves on to the next line
        }
        
        $line =~ s/\R//g;  #removes any linebreaks from the data line being read (including the one at the end)
        $progress ++;  #Increments the progress (lines processed) count
        my @line = split(/\t/, $line);  #splits the data line just like the header
        if (scalar(@line) != scalar(@header)) {  #tests to see if there as many data entries as there are headers
            print "\nMismatched header and data on line $progress.\n";  #Prints the error to the console
            print LOGFILE "\nMismatched header and data on line $progress.\n"; #Prints the error to the logfile
            next LINE;  #skips the line if there are not
        }
        
        my $index = 0; #resets the index of which column within the row is being moved to the hash
        my %linehash = ();  #creates a hash with the header entry as key and data as value
        
        foreach my $headervalue(@header){  #goes through each header entry
            $linehash{$header[$index] } = $line[$index];  #writes the header entry as key and data entry as value in the hash
            $index++;  #increments the index so that the next iteration will look at the next column in each line
        }
        
        my $millimorgan = int($linehash{position}/100000);  #calculates which 100000 base segment is being checked to select the correct reference file
        
        my $reference = "subvcfs\/".$linehash{chromosome}."c".$millimorgan.".subvcf"; #generates the filename to open for a snp reference
        
        my @refline;  # initializes the array variable containing all the reference line data
        #my $refpos = 0;
        my $allpop5 = "False";  # initializes the variable to false
        my $onepop5 = "False";  # initializes the variable to false
        my $mut = "False";  # initializes the variable to false
        my $rsid = "None";  # initializes the variable to none
        my $obsallele = $linehash{sampleAlleles};
        #my $key = 0;
        $obsallele =~ s/\///g; # removes any slashes from the observed allele field
        my @obsallele = split(//, $obsallele);  #splits the observed allele into an array of characters

        unless (-e "$reference"){  #checks if the SNP reference file exists
            $lastmillimorgan = $millimorgan;  #if it does not, assumes there are not SNPs reproted in the region 
            print OUTPUT $line."\t".$rsid."\t".$allpop5."\t".$onepop5."\t".$mut."\n"; #writes the output to the file
            print "\rProcessed $progress lines\."; #updates the progress report
            next LINE; #goes on to next line
        }
        
        unless ($millimorgan == $lastmillimorgan){  #checks to see if the current variant to be annotated is in the same block of 100000 bases as the last
            open(SNPFILE, $reference)  or die "\nMissing subvcf file $reference\.";  #if not, opens the appropriate reference file
            local $/ = undef;  #tells the computer to ignore any end of line characters and read the whole file as a line
            $snplibrary = <SNPFILE>; #actually tells it to read the file into a single variable
            close SNPFILE;  #closes the snp reference file
        }
        
        my $regex = "\(".$linehash{chromosome}."\t".$linehash{position}."\t.*?\)\n";  #creates a string of text to search the variable containing the reference by chromosome and position
        
        unless ($snplibrary =~ /$regex/i){  #checks to see if a matching chromosome and position can be found in the reference file
            $lastmillimorgan = $millimorgan;  #if not (because the locus is not listed in dbSNP), remembers the current block of 100,000 bases for the next line
            print OUTPUT $line."\t".$rsid."\t".$allpop5."\t".$onepop5."\t".$mut."\n";  #prints the annotation to the output file
            print "\rProcessed $progress lines\.";  #updates the progress line
            next LINE; #moves on to the next variant
        }
        
        
            my $refline = $1;  #captures the entire line of what was found when reading the reference file containing the appropriate locus (this is done using the parentheses in $regex)
            @refline = split (/\t/, $refline);   #splits the reference line on each tab into an array 
            $rsid = $refline[2];  #perl arrays are indexed (their first value) to zero, so this saves the third value in the array as the rs number
            unless ($linehash{referenceBase} eq $refline[3]) {  #checks to see if the reference bases in the seattleseq output and dbSNP are the same
                print "\nReference bases at $refline[0]\:$refline[1] do not match\, check reference genome versions in both files\.\n ";  #if not, prints an error message
                print LOGFILE "\nReference bases at $refline[0]\:$refline[1] do not match\, check reference genome versions in both files\.\n ";  #then prints the annotated line to the output file
                $allpop5 = "ref base mismatch";  #annotates that ref bases don't match
                $onepop5 = "ref base mismatch";  #same as above
                $mut = "ref base mismatch";  #same as above
                $lastmillimorgan = $millimorgan;  #remembers the current 100000 base block
                print OUTPUT $line."\t".$rsid."\t".$allpop5."\t".$onepop5."\t".$mut."\n";  #prints the annotated line to the output file
                print "\rProcessed $progress lines\.";  #updates the progress report
                next LINE;  #goes to the next line
                }
                
            foreach my $allele(@obsallele){  #checks each of the observed alleles
                if ($allele eq $refline[4]) {  #asks if the observed allele matches with the dbSNP variant
                  if ($refline =~ /\;MUT/i) {  #if so, checks if the line is annotated as a mutation
                    $mut = "True";  #and returns a true value for $mut if it is
                  }
                  if ($refline =~ /\;G5A/i) {  #same as above for G5A (if the variant is over 5% frequency in all populations)
                    $allpop5 = "True";
                  }
                  if ($refline =~ /\;G5\;/i) {  #same as above for G5 (if the variant is over 5% frequency in at least 1 population)
                    $onepop5 = "True";
                  }
                else{  #if not...
                    if ($allele ne $refline[3]) {  #if the allele is not the same as the reference
                        $allpop5 = "False";  #annotates it as something unknown in dbSNP
                        $onepop5 = "False";  #annotates it as something unknown in dbSNP
                        $mut = "False";  #annotates it as something unknown in dbSNP
                        $lastmillimorgan = $millimorgan;  #remembers the current block of 100000 bases for the next line
                        print OUTPUT $line."\t".$rsid."\t".$allpop5."\t".$onepop5."\t".$mut."\n";  #prints the output
                        print "\rProcessed $progress lines\.";  #updates progress
                        next LINE;  #moves on to the next line
                        }
                        
                      }
                                         
                    }
                }
                
            
            
            
        
        $lastmillimorgan = $millimorgan;  #remembers the current block of 100000 for the next line
        print OUTPUT $line."\t".$rsid."\t".$allpop5."\t".$onepop5."\t".$mut."\n";  #prints the annotation to the output file
        print "\rProcessed $progress lines\.";  #updates the progress counter
    }
    close OUTPUT;  #closes the output file
    close LOGFILE;  #closes the log file
    print "\nDone\!\n";  #prints done in the console
}

sub checkopts{
    my $opts = shift;  #dereferences the hash containing the options
    
    my $file = $opts->{"f"}; #puts the value in options under key F into a variable called file
    
    unless(defined($file) and (-e $file)){  #unless the file entered exists...
        print "Input file not found or not defined in commandline arguments.\n";
        return 0;  #this function will return a value of 0, which is false and signals bad options
    }
}

sub usage{  #This subroutine prints directions
    print "This program will append additional dbSNP data on a SeattleSeq output.\nSample commandline\:\nperl dbSNPerator\.pl \-f file\.vcf\.\n";
    die "";
}

sub checkintegrity{  #This subroutine checks make sure files are not missing from the subvcf directory
    unless(opendir(SUBVCFDIR,"subvcfs")){  #If the subvcf directory cannot be opened
        die "\nError opening subvcf directory\.\n";  #program quits
    }
    my $filelist;  #initiates a variable called filelist
    my @files = readdir("SUBVCFDIR");  #creates an array of all the file names in the subvcf directory
    foreach my $file(@files){  #for all of the filenames in the subvcf directory
        $filelist = $filelist.$file;  #add the filename to the growing string of filenames
    }
    closedir(SUBVCFDIR);  #closes the subvcf directory
    my $dirdigest = md5($filelist);  #creates an md5 digest of the filename list string
    open (CHECK, 'subvcfs/0test.subvcf') or die ("\nUnable to create integrity check file\.\n");
    my $stored = <CHECK>;
    close CHECK;
    chomp ($stored);
    return $dirdigest eq $stored;
    }


main();
