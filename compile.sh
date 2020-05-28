echo Compling SSumM java source codes...
rm -rf class
mkdir class

javac -cp ./fastutil-8.2.3.jar -d class $(find ./src -name *.java)

echo Make SSumM jar archive...
cd class
jar cf SSumM.jar ./
rm ../SSumM.jar
mv SSumM.jar ../
cd ..
rm -rf class
cd ..

echo Done.
