#!/usr/bin/perl
open(FH,">> DS.dat");
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1600 200`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1840 184`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 2040 204`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1600 200`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 1840 184`;
    print FH $tx;
}
for ($mn = 4000; $mn <= 40000; $mn = $mn + 2000) {
    $tx = `./DS $mn $mn 2040 204`;
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
