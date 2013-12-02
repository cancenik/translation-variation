## CC 
## May 20 2013

# This script runs windowFilter sequentialls on all files in a directory

for f in *
do
    windowFilter.pl -i $f -o WINDOW_OUT -m 1 -c 6
done
