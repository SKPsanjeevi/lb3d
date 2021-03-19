# @ job_name = JUQUEEN_job
# @ comment = "JUQUEEN job"
# @ error = $(job_name).$(jobid).err
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL
# @ wall_clock_limit = 0:30:00
# @ notification = always
# @ notify_user = first.lastname@your.site.edu
# @ job_type = bluegene
# @ bg_size = 128
# @ queue

# Use 'llsubmit <jobfile>' to submit.

runjob --ranks-per-node 16 : ./lbe -f input-file
