wget https://adlibitum.oats.inaf.it/monaco/etc/perAbInf.tgz
tar zxf perAbInf.tgz
echo path of the folder to be created
read path
mkdir -p $path
cd data
for DIR in *
do
 cp -R $DIR $path
done

