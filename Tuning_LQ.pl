#!/usr/bin/perl
open(FH,">> LQ.dat");
foreach $mn ( 8000, 15000, 15500, 16000, 16500, 17000 ) {
      $tx = `./LQ $mn $mn 2`;
      print FH $tx;
}
close(FH);
