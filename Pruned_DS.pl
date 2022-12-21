#!/usr/bin/perl
open(FH,">> DS_201602.dat");
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 160 80`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 320 80`;
    print FH $tx;
}
#for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
#    $tx = `./DS $mn $mn 480 120`;
#    print FH $tx;
#}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 720 120`;
    print FH $tx;
}
#for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
#    $tx = `./DS $mn $mn 960 120`;
#    print FH $tx;
#}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1200 200`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1440 240`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1680 240`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1920 240`;
    print FH $tx;
}
close(FH);
$sendmail = '/usr/sbin/sendmail';
$from = 'stomo@yamanashi.ac.jp';
$to = 'stomo@yamanashi.ac.jp';
$subject = 'DS Finish';
$msg = 'DS Finish';
open(SDML,"| $sendmail -t -i") || die 'sendmail error';
print SDML "From: $from\n";
print SDML "To: $to\n";
print SDML "Subject: $subject\n";
print SDML "$msg";
close(SDML);
