rm SSumM.tar.gz
rm -rf SSumM
mkdir SSumM
cp -R ./{compile.sh,run.sh,package.sh,ego_facebook.txt,fastutil-8.2.3.jar,LICENSE,src,Makefile,README.txt,user_guide.pdf,output} ./SSumM
tar cvzf SSumM.tar.gz --exclude='._*' --exclude="*.iml" ./SSumM
rm -rf SSumM
echo done.
