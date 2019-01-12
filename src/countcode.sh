echo "MIRA"
find errorhandling io lec mira stdinc util -name "*.C" >ttt 
find errorhandling io lec mira stdinc util -name "*.H" >>ttt 
echo "Files: " `grep -c ^ ttt`
find errorhandling io lec mira progs stdinc util -name "*.C" -exec cat {} >ttt \;
find errorhandling io lec mira progs stdinc util -name "*.H" -exec cat {} >>ttt \;
grep -v '^[[:space:]]*[{}]*[[:space:]]*$' ttt >ttt2
grep -v '^[[:space:]]*//.*$' ttt2 > ttt
echo "Classes: " `grep ^class ttt`
echo "Classes: " `grep -c ^class ttt`
echo "LOC Count: " `wc ttt`
echo

echo "EdIt"
find EdIt caf examine -name "*.C" >ttt
find EdIt caf examine -name "*.H" >>ttt
echo "Files: " `grep -c ^ ttt`
find EdIt caf examine -name "*.C" -exec cat {} >ttt \;
find EdIt caf examine -name "*.H" -exec cat {} >>ttt \;
grep -v '^[[:space:]]*[{}]*[[:space:]]*$' ttt >ttt2
grep -v '^[[:space:]]*//.*$' ttt2 > ttt
echo "Classes: " `grep -c ^class ttt`
echo "LOC Count:" `wc ttt`
rm ttt ttt2
