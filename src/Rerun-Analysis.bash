echo "Running the whole analysis can take days."

bash Pre1-Trimming-Reads.bash
bash Pre2-Mapping-Reads-Subread.bash
bash Pre3-Mapping-Reads-STAR.bash
bash Tab1-Mapped-Reads.bash
bash Tab2-Coef-with-Truth.bash
bash SupplTab1-Trimmed-Bases.bash
bash SupplTab2-Adapter-Clipped.bash
