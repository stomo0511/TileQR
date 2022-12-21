#!/usr/bin/perl
open(FH,">> KP.dat");
foreach $mn ( 2000, 4000, 6000, 8000, 10000, 12000, 14000 ) {
#  foreach $bs ( 50, 100, 200, 300, 400, 500, 600 ) {
  foreach $bs ( 300, 400 ) {
#    foreach $ib ( 10, 25, 50, 100, 200, 300 ) {
    foreach $ib ( 300 ) {
      if ($bs >= $ib) {
	  for ( $i=0; $i<5; $i++) {
	      $tx = `./KL $mn $mn $bs $ib`;
	      print FH $tx;
	  }
      }
    }
  }
}
close(FH);
