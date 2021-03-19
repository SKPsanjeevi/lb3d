#/bin/bash
# This script searches all .F90 and .h files for #ifdef <FLAG> statements and reports all <FLAG>s found.
#
# ===== Copy this file to the code directory to execute it! =====
#

# Write two temp files - one with just <FLAG>, the other with <FLAG> <TOTAL COUNT>
ls *.F90 *.h | xargs grep '^#e\?l\?ifn\?\(def\)\?' | sed 's/.*#e\?l\?ifn\?\(def\)\? \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq > find_flags.tmp.total
ls *.F90 *.h | xargs grep '^#e\?l\?ifn\?\(def\)\?' | sed 's/.*#e\?l\?ifn\?\(def\)\? \(.*\)/\2/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > find_flags.tmp.total.cnt

# Mostly the same as above, except we only look inside the lbe_detect_flags() subroutine
cat lbe_io.F90 | awk '/subroutine lbe_detect_flags()/,/end subroutine lbe_detect_flags/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq > find_flags.tmp.lbe_detect_flags
cat lbe_io.F90 | awk '/subroutine lbe_detect_flags()/,/end subroutine lbe_detect_flags/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > find_flags.tmp.lbe_detect_flags.cnt

# Mostly the same as above, except we only look inside the lbe_more_metadata_phdf5() subroutine
cat lbe_io_hdf5.F90 | awk '/subroutine lbe_more_metadata_phdf5()/,/end subroutine lbe_more_metadata_phdf5/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq > find_flags.tmp.lbe_more_metadata_phdf5
cat lbe_io_hdf5.F90 | awk '/subroutine lbe_more_metadata_phdf5()/,/end subroutine lbe_more_metadata_phdf5/' | grep '^#ifdef' | sed 's/.*#ifdef \(.*\)/\1/' | sed 's/[ \t]*$//' | sort | uniq -c | awk '{print$2, $1}' > find_flags.tmp.lbe_more_metadata_phdf5.cnt

# Add <FLAG> 0 entries for the missing flags
cat find_flags.tmp.total find_flags.tmp.lbe_detect_flags | sort | uniq -u | sed 's/\(.*\)/\1 0/' >> find_flags.tmp.lbe_detect_flags.cnt
cat find_flags.tmp.total find_flags.tmp.lbe_more_metadata_phdf5 | sort | uniq -u | sed 's/\(.*\)/\1 0/' >> find_flags.tmp.lbe_more_metadata_phdf5.cnt

# Sort the files so the entries exactly match the order of the main file
sort find_flags.tmp.lbe_detect_flags.cnt > find_flags.tmp.lbe_detect_flags.cnttotal
sort find_flags.tmp.lbe_more_metadata_phdf5.cnt > find_flags.tmp.lbe_more_metadata_phdf5.cnttotal

# Just some text written to the report
echo "Current compiler flags:" > find_flags_report.txt
echo "NAME #total #stdout #hdf5" >> find_flags_report.txt
echo " " >> find_flags_report.txt

# Paste the lines of all files, line by line, and retain only <FLAG> <TOTAL COUNT> <IO> <IO_HDF5>
paste find_flags.tmp.total.cnt find_flags.tmp.lbe_detect_flags.cnttotal find_flags.tmp.lbe_more_metadata_phdf5.cnttotal | awk '{print$1, $2, $4, $6}' >> find_flags_report.tmp

# Cleanup
rm find_flags.tmp.*

# Append data to final report
cat find_flags_report.tmp >> find_flags_report.txt
echo " " >> find_flags_report.txt
echo "Checking for errors: " >> find_flags_report.txt

echo "Writing report to find_flags_report.txt"
echo "Analyzing..."

# Line by line, check if the data is what we want
cat find_flags_report.tmp | while read line; do
  varn=`echo $line | awk '{print$1}'`
  vart=`echo $line | awk '{print$2}'`
  ldf=`echo $line | awk '{print$3}'`
  lmmp=`echo $line | awk '{print$4}'`
  # Total number of entries outside the output routines
  vart=$(( $vart - $ldf - $lmmp ))
  # We need exactly one entry for <IO>
  if [ $ldf -ne 1 ]; then
    echo "WARNING: Flag $varn not reported exactly once in lbe_detect_flags()." | tee -a find_flags_report.txt
  fi
  # We need exactly one entry for <IO_HDF5>
  if [ $lmmp -ne 1 ]; then
    echo "WARNING: Flag $varn not reported exactly once in lbe_more_metadata_phdf5()." | tee -a find_flags_report.txt
  fi
  # We need at least one entry outside those two output routines
  if [ $vart -lt 1 ]; then
   echo "WARNING: Flag $varn is reported but unused." | tee -a find_flags_report.txt
  fi
done

# More text output to the report
echo " " >> find_flags_report.txt
echo "Done!" >> find_flags_report.txt

# Cleanup
rm find_flags_report.tmp

echo "Done!"

