#!/bin/csh
##########################################################
# Junhong Lee / Shell for Windprofiler data / 160908
##########################################################
# authkey ; junhong : efbc706c7892b012a0363e3c3b20548c073e06a8e9a0c0d554003653725538b1c0f3263d6b910ee3e176fdd61a3121ed6184c373ed823064a51eaf6120893447

set dir = "/data1/junhong/WRF/pblh/AWS/data"
set key = "efbc706c7892b012a0363e3c3b20548c073e06a8e9a0c0d554003653725538b1c0f3263d6b910ee3e176fdd61a3121ed6184c373ed823064a51eaf6120893447"

cd $dir
cp ../url_test2 .

#-------------- collection period & Station number --------------------
#             ( YYYY MM DD HH ) # 10 MIN interval data... !! 
set start   = ( 2016 04 28 18 )
set end     = ( 2016 05 01 00 )

#----------------------------------------------------------------------
set startday = $start[1]$start[2]$start[3]$start[4]
set endday  = $end[1]$end[2]$end[3]$end[4]

set lyear   = $start[1]
set lmonth  = $start[2]
set lday    = $start[3]
set lhour   = $start[4]
set ndate = $startday

while ( $ndate <= $endday )

#     set ryear   = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%Y'`
#     set rmonth  = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%m'`
#     set rday    = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%d'`
#     set rhour   = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%H'`
#     set rdate = $ryear$rmonth$rday$rhour

#    gcc -c $dir/url_io.c $dir/url_fort.c
#    g77 -o $dir/url_test2 $dir/url_test2.f $dir/url_io.c $dir/url_fort.o -lcurl

#    ./url_test2 "http://203.247.66.28/url/cgi-bin/url/nph-aws2_min_raw?tm1=${sdate}&tm2=${ndate}&stn=${station}&disp=0&help=1&authKey=a8ecd9c7c12572c4596ecbd5b7e13b5ed15a2cf8a5492bba5642f293ee63758f099dec79b3d7158174c22994f473144c6f74b61c7bb47db6302247dba517d202">>'AWS2_'$station'_'${startday}'_'${endday}.dat
#./url_test2 "http://203.247.66.28/url/kma_wpf.php?tm=${ndate}&stn=${station}&mode=L&help=0&authKey=${key}" > UPPER_WP_ALLST_ULC_${station}_${ndate}.raw.data
#./url_test2 "http://203.247.66.28/url/kma_wpf.php?tm=${date}&stn=0&mode=L&help=0&authKey=6b89a8e82e7dc3a444ec0d0b0dbbd21ef86bc5c2151e15f9a5629a3ca147bb56cda47ea0a5305014474e1404bd68f65128fb2b80ce121bc85e4909c0b5c175a9"


./url_test2 "http://203.247.66.28/url/awsh.php?var=WD&tm=${ndate}00&help=0&authKey=${key}" > ${ndate}_hourl_wind.dat

#echo "doing" $ndate UPPER_WP_ALLST_ULC_${station}_${ndate}.raw.data
lo ${ndate}_hourl_wind.dat

     set nyear   = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%Y'`
     set nmonth  = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%m'`
     set nday    = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%d'`
     set nhour   = `date -d "${lyear}-${lmonth}-${lday} ${lhour}:00:00 1 hours" '+%H'`
     set lyear   = $nyear
     set lyear   = $nyear
     set lmonth  = $nmonth
     set lday    = $nday
     set lhour   = $nhour
     set ndate = $nyear$nmonth$nday$nhour

end

echo "data from" $startday "to" $endday "!"

exit
