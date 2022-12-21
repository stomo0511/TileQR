#!/usr/bin/perl
open(FH,">> RT_80_40.dat");
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./RT $mn $mn 80 40`;
    print FH $tx;
}
close(FH);
