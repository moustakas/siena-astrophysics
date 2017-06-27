#!/bin/bash
# Joey Rowley, Kenneth Tousignant, John Moustakas, & John Cummings
# 2014 Jun

wwwdir='/var/www/siena_weather_cam/'
allskydir='/data/allskycam/'

today=$(date +%Y%m%d)
yesterday=$(date +%Y%m%d -d "yesterday")
latestjpg_www='allskycam-latest.jpg'
latestpanojpg_www='allskycam-panoramic.jpg'
movieformat='.mov'

suffix='????.jpg'
panosuffix='????-panoramic.jpg'

#makes directory for the current day
echo "Making directory $allskydir$today"
mkdir -p $allskydir$today

#moves the files to correct directory
echo "Moving files $allskydir$today$suffix to $allskydir$today"
mv -f $allskydir$today$suffix $allskydir$today 2> /dev/null
mv -f $allskydir$today$panosuffix $allskydir$today 2> /dev/null

#updates permissions
echo "Updating permissions."
chmod 775 $allskydir$today
chmod 664 $allskydir$today/$today$suffix
chmod 664 $allskydir$today/$today$panosuffix

chgrp -R weather $allskydir$today

#creates solft link that links to the website directory
latestjpg=`ls -tr $allskydir$today/$today$suffix | tail -1`
latestpanojpg=`ls -tr $allskydir$today/$today$panosuffix | tail -1`
echo "ln -sf $latestjpg $wwwdir$latestjpg_www"
echo "ln -sf $latestpanojpg $wwwdir$latestpanojpg_www"
ln -sf $latestjpg $wwwdir$latestjpg_www
ln -sf $latestpanojpg $wwwdir$latestpanojpg_www

# --------------------------------------------------
# make the last-hour movie

# get all of today's and yesterday's images
yesterfiles=`ls $allskydir$yesterday/${yesterday}????.jpg 2> /dev/null`
todayfiles=`ls $allskydir$today/${today}????.jpg 2> /dev/null`
allfiles=(${yesterfiles[@]} ${todayfiles[@]})
#echo ${allfiles[@]} # here

# get the year-month-day-hour-minute string corresponding to one hour ago
#onehourago='201406110033'
onehourago=$(date +%Y%m%d%H%M -d "1 hour ago")
#echo $onehourago # here

# loop on each image and find those files whose time stamp indicates
# it was obtained more than ONEHOURAGO
for file in ${allfiles[@]}; do
    filehour=$(basename $file ".jpg")
    if [ $filehour \> $onehourago ];
    then
	lasthourfiles=(${lasthourfiles[@]} $file)
#	echo $thishour $onehourago;
    fi;
done;
#echo ${lasthourfiles[@]} # here

# build the movie by copying the relevant images to a temporary
# directory before calling ffmpeg
tmpdir=`mktemp -d`
#echo "Making temporary directory $tmpdir" # here

let i=0
for file in ${lasthourfiles[@]}; do
#   echo "cp $file $tmpdir/$(printf %05d.jpg ${i})"
    cp $file $tmpdir/$(printf %05d.jpg ${i})
 	let i++
done;

lasthourmovie=$allskydir/allskycam-lasthour$movieformat
ffmpeg -r 3 -i $tmpdir/%05d.jpg -y $lasthourmovie
#rm -rf $tmpdir # here

# update permissions and the link on the website to the latest movie
chmod 664 $lasthourmovie
chgrp weather $lasthourmovie
ln -sf $lasthourmovie ${wwwdir}/allskycam-lasthour$movieformat

# --------------------------------------------------
# make the all-day movie, where the "day" starts at 6pm

# get all the images between 6pm yesterday and 5:5pm today
yesterfiles1=`ls $allskydir$yesterday/${yesterday}1[8-9]??.jpg 2> /dev/null`
yesterfiles2=`ls $allskydir$yesterday/${yesterday}2[0-3]??.jpg 2> /dev/null`
yesterfiles=(${yesterfiles1[@]} ${yesterfiles2[@]})
#echo ${yesterfiles[@]}

todayfiles1=`ls $allskydir$today/${today}0???.jpg 2> /dev/null`
todayfiles2=`ls $allskydir$today/${today}1[0-7]??.jpg 2> /dev/null`
todayfiles=(${todayfiles1[@]} ${todayfiles2[@]})
#echo ${todayfiles[@]}

# combine the file lists
alldayfiles=(${yesterfiles[@]} ${todayfiles[@]})
#echo ${alldayfiles[@]}

# build the movie by copying the relevant images to a temporary
# directory before calling ffmpeg
tmpdir=`mktemp -d`
#echo "Making temporary directory $tmpdir"

let i=0
for file in ${alldayfiles[@]}; do
#   echo "cp $file $tmpdir/$(printf %05d.jpg ${i})"
    cp $file $tmpdir/$(printf %05d.jpg ${i})
    let i++
done;

alldaymovie=$allskydir/${today}/allskycam-${today}$movieformat
ffmpeg -v 0 -r 10 -i $tmpdir/%05d.jpg -y $alldaymovie
rm -rf $tmpdir

# update permissions and the link on the website to the latest movie
chmod 664 $alldaymovie
chgrp weather $alldaymovie
ln -sf $alldaymovie ${wwwdir}/allskycam-today$movieformat
