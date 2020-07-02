
if [[ ! ( ` python -V 2>&1 ` =~ "Python 3" ) ]]
then
	echo "You need to have Python 3 in PATH and its name should be 'python' (not python3)."
	exit
fi
echo "Running the whole analysis can take days. "
read -p "Press ENTER to continue, or ^C to stop." del4

bash Pre1-Trimming-Reads.bash
bash Pre2-Mapping-Reads-Subread.bash
bash Pre3-Mapping-Reads-STAR.bash
bash Tab1-Mapped-Reads.bash
bash Tab2-Coef-with-Truth.bash
bash SupplTab1-Trimmed-Bases.bash
bash SupplTab2-Adapter-Clipped.bash
