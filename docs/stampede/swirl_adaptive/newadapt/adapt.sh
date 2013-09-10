#! /bin/sh

cat adapt.txt | \
awk '{ printf ("%.1f %.1f %.1f\n", $1, $2, $1 / $2 / 16 * 100); }' \
> adapt.tab
