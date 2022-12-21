#!/usr/bin/perl
open(FH,">> error.dat");
$mn = 30720;
$tx = `./DS $mn $mn 80 16`;
print FH $tx;
$tx = `./DS $mn $mn 160 32`;
print FH $tx;
$tx = `./DS $mn $mn 320 64`;
print FH $tx;
$tx = `./DS $mn $mn 640 128`;
print FH $tx;
$tx = `./DS $mn $mn 1280 256`;
print FH $tx;
close(FH);
