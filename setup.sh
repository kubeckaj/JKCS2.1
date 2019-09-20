# USE
# ./setup.sh
# ./setup.sh -r //rewrite ~/.JKCSusersetup.txt

# START
echo "-----------------------"
echo 'Hi, I will just test what is written in your ~/.bashrc file'

# PATH DECLARATION
path=$PWD
path1="${path}/JKCSx"
path2="${path}/TOOLS/MANIPULATE"

# DOES ~/.bashrc exist?
if [ ! -e ~/.bashrc ]
then
  touch ~/.bashrc
fi

function writetobashrc {
  command="$1"
  test=`grep -c "$command" ~/.bashrc`
  if [ $test -ne 0 ]
  then 
    sed  "\#$command#d" ~/.bashrc > ~/.bashrc_help
    mv ~/.bashrc_help ~/.bashrc
  fi
  echo "$command" >> ~/.bashrc
}

command1="export PATH=$path1:\$PATH"
writetobashrc "$command1"

command2="export PATH=$path2:\$PATH"
writetobashrc "$command2"

echo "Write followinng command: "
echo "          source ~/.bashrc"
echo "-----------------------"

if [ ! -e ~/.JKCSusersetup.txt ] || [ "$1" == "-r" ]
then
  if [ -e ~/.JKCSusersetup.txt ]; then cp ~/.JKCSusersetup.txt ~/.oldJKCSusersetup.txt; fi
  cp .JKCSusersetup.txt ~/.JKCSusersetup.txt
  echo "Please, change all required paths in file ~/.JKCSusersetup.txt"
else
  echo "File ~/.JKCSusersetup.txt already exist, but check, if all paths in this file are correct. (or use -r argument to rewrite)"
fi
echo "-----------------------"
