# @Author: Daniel
# @Date:   2017-05-12 10:32:53
# @Last Modified by:   Daniel
# @Last Modified time: 2017-10-27 11:44:32


source ~/.bash_profile
wd=`pwd`
echo "pwd: " `pwd`
echo "config: $wd/$1"
cmd='Arguments='"$wd/$1"' --jobIndex=$(Process)'
condor_submit grid/config.condor.sh -append "$cmd"

