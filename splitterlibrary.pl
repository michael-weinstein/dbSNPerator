#!/usr/bin/perl -w

#Cohn lab script for splitting a large VCF into smaller reference libraries.  Required for the dbSNPerator script to function.
#Copyright 2014 Michael Weinstein, Daniel H. Cohn laboratory, UCLA.

use strict;
use warnings;
use Getopt::Std;
use Digest::MD5 qw(md5 md5_hex md5_base64);
$|=1;

sub main{
    my %opts;  #creates a hash table for options
    getopts('s:', \%opts);  #actually recalls commandline options and enters them into the hash
    unless(checkopts(\%opts)){  #unless the entered options were valid
        usage();  #runs the usage subroutine
    }
    
    my $file = $opts{s};  #looks at the options to find the filename
    
    print "\nCreating subvcf directory\.\n";  #notifies user of subdirectory creation
    if (-d "subvcfs") {die "Subvcf directory already exists\!\n";} #quits if there is already a subdirectory by that name
    mkdir('subvcfs'); #actually makes the subdirectory
    
    open (FILETOSPLIT, $file) or die ("Error opening file $file.\n");  #opens the large snp reference VCF
    
    my $progress = 0;  #initializes a variable to monitor progress
    my $currentchr = 0;  #initializes a variable value for the current chromosome
    my $millimorgan = 0;  #initializes a variable value for the current block of 100000 bases
    
    print "Starting the split \(this may take a while\)\.\n\n";  #notifies user that the split is starting
    
    READLINE: while (my $line = <FILETOSPLIT>){  #reads a line from the SNP reference VCF
        print "\rProcessed $progress lines."; #Gives an output of the progress
        if ($line =~ /^#/){ #if the line starts with a # symbol
            $progress ++; #increment the progress counter
            next READLINE;  #go to the next line in the reference
            }
        if ($line =~ /^\d+?|x|y|mt*?\t/i) {  #if the line starts off with a number, the letters X, Y, or MT
            my @line = split(/\t/, $line);  #splits the line into an array of values
            $currentchr = $line[0];  #writes the current chromosome to the appropriate variable
            $millimorgan = int($line[1]/100000);  #writes the current block of 100000 bases to the appropriate variable
            my $chrfile = "subvcfs\/".$line[0]."c".$millimorgan.".subvcf"; #generates a name for the appropriate library file to edit, looks like [chromosome]c[block of 100000 bases]
            open (FILEOUT, '>>'.$chrfile) or die ("\nUnable to open output file $chrfile.\n"); #opens that file
            print FILEOUT $line; #writes the reference snp VCF line to the smaller library file
            close FILEOUT;  #closes the library file
            }
 
        
        else { #if the line doesn't fit either of those criteria
            my $badline = $progress + 1;   #notes the line number where the error occurred
            print "\nError processing line $badline\. Does not appear to start with a chromosome number.\nReadout\: $line\nNote that a blank line can often be found at the end of a file\.\n"  #notifies user of a possible bad line
            }
                   
                   
        $progress ++;  #increments the progress counter
        next READLINE;  #moves on to the next line in the reference VCF
    }
    close FILETOSPLIT;  #after all is done, closes the reference VCF
   
}
    
sub checkopts{ #subroutine to check options
    my $opts = shift;  #dereferences the hash of options
    
    my $file = $opts->{"s"}; #captures the filename entered in commandline options
    
    unless(defined($file) and (-e $file)){ #unless the filename was entered and such a file exists
        return 0; #return a value of 0 (false)
    }
}

sub usage{  #prints the directions for usage and then quits
    print "This program will split a large VCF into sub VCFs by chromosome for easier comparison.\nSample commandline\:\nperl splitter\.pl \-s file\.vcf\.\n";
    die "";
}

sub integrity{  #integrity checking subroutine to generate the md5 hash for confirmation in dbsnperator
        unless(opendir(SUBVCFDIR,"subvcfs")){ #if the directory can't be opened,
            die "\nError opening subvcf directory\.\n";  #quits the program and displays an error message
    }
    open (DIGESTOUT, '>>subvcfs/0test.subvcf') or die ("\nUnable to create integrity check file\.\n"); #creates a file for saving the md5 hash value
    my $filelist; #initializes a string for saving the list of file names
    my @files = readdir("SUBVCFDIR"); #creates an array of filenames from the library directory
    foreach my $file(@files){  #for all the values in the filename array
        $filelist = $filelist.$file; #add it to a (growing) string containing all the filenames
    }
    closedir(SUBVCFDIR); #closes the directory
    my $dirdigest = md5($filelist); #creates an md5 digest of the long string of filenames
    
    print DIGESTOUT "$dirdigest"; #prints the value of the digest to the file where it is being stored
    close DIGESTOUT; #closes the file where the digest is stored
}

main();  #actually runs the main program
integrity(); #actually runs the md5 hash generating integrity subroutine


