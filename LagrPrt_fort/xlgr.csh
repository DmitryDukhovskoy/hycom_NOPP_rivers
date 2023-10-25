#!/bin/csh -vx
make clean
make
wait
if (! $status ) then
  ./lagr.x
endif

exit 0


