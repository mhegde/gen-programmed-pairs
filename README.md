
This script helps make sgRNA combinations for programmed pairs of genes.
<p><b>Inputs</b>
<ol> 
<li>S.pyogenes guide designs: .csv file with S.pyogenes guide designs for all genes</li>
<li>S.aureus guide designs: .csv file with S.pyogenes guide designs for all genes</li>
<li>S.pyogenes control guides: .txt file with list of S.pyogenes controls. Column name should be 'sgRNA Seq'</li>
<li>S.aureus control guides: .txt file with list of S.aureus controls. Column name should be 'sgRNA Seq'</li>
<li>Gene pairs: .txt file with gene pairs. Refer to example file provided for format</li>
<li>Overlap check: Number of nucleotides overlap between sgRNA pairs to be avoided. Default:12</li>
<li>Outputfile: Name of outputfile</li>
