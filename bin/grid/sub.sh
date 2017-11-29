# @Author: Daniel
# @Date:   2017-05-12 10:32:53
# @Last Modified by:   Daniel
# @Last Modified time: 2017-10-27 11:44:32


source ~/.bash_profile
wd=`pwd`
echo "pwd: " `pwd`
echo "config: $wd/$1"

cmd1="InitialDir=$wd"
cmd2="Executable=$wd/pairAna.app"
cmd3='Arguments='"$wd/$1"' --jobIndex=$(Process)'

echo "==============="
echo $cmd1
echo $cmd2
echo $cmd3
echo "==============="

condor_submit grid/config.condor.sh -append "$cmd1" -append "$cmd2" -append "$cmd3" --queue $2

