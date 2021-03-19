# @ job_name = Single
# @ comment = "BGP runny hopphopp"
# @ environment = COPY_ALL;
# @ wall_clock_limit = 0:25:00
# @ notification = always
# @ notify_user = jens@icp.uni-stuttgart.de
# @ job_type = bluegene
# number of processors:
# @ bg_size = 32
# @ error = $(job_name).$(bg_size).$(jobid).err
# @ output = $(job_name).$(bg_size).$(jobid).out
# @ queue
/usr/local/bin/mpirun -exe `/bin/pwd`/lbe \
 -verbose 1 -args " -f input-test2"


