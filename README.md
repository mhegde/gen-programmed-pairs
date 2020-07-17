
This script helps make sgRNA combinations for programmed pairs of genes.
<p><b>Inputs</b>
<ol> 
<li>S.pyogenes guide designs: .csv file with S.pyogenes guide designs for all genes</li>
<li>S.aureus guide designs: .csv file with S.pyogenes guide designs for all genes</li>
<li>S.pyogenes no-site guides: .txt file with list of S.pyogenes no-site controls. Column name should be 'sgRNA Sequence' and 'Target Gene Symbol'</li>
<li>S.pyogenes one-non-gene-site guides: .txt file with list of S.pyogenes one_non-gene_site controls. Column name should be 'sgRNA Sequence' and 'Target Gene Symbol'</li>
<li>S.aureus no-site guides: .txt file with list of S.aureus no-site controls. Column name should be 'sgRNA Sequence' and 'Target Gene Symbol'</li>
<li>S.aureus one-non-gene-site guides: .txt file with list of S.aureus one_non-gene_site controls. Column name should be 'sgRNA Sequence' and 'Target Gene Symbol'</li>
<li>Gene pairs: .txt file with gene pairs. Refer to example file provided for format</li>
<li>Overlap check: Number of nucleotides overlap between sgRNA pairs to be avoided. Default:12</li>
<li>Number of control guides with targeting guides: Number of control guides to be paired with targeting guides. Same number of no-site and intergenic sites. Default:5</li>
<li>Number of control guides with control guides: Number of control guide combinations. Same number of no-site and intergenic sites. Default:50</li>
<li>Outputfile: Name of outputfile</li>
