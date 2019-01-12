
FROM=strainname
TO=strain_name

for f in *xml; do
  echo $f
  sed -e "s/$FROM/$TO/g" $f >bla
  mv bla $f
done
