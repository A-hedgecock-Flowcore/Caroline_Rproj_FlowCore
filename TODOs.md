[x] add cleaning algorithm before gating steps (PeacoQC)?
[x] add a quarto doc that walks thru manually editing gates if necessary
[] fix compensAID or just stick with NxN
[x] DR and clustering script
 [x]   - separate after QC script is verified
[x] pdf version of QC report
[] annotating clusters dynamically?
[] should controls be processes the same? (LiveDead only and PBMC + U251 control)
[] gating file
[] transformation of certain poorly separated markers (i.e. CD4)
[] update README

[] For the cleanup gating, do you think we could change the order to Helix+-> CD45+ -> singlets -> Live/Dead?
That is more how I do it in the manual analysis to make the live dead gating less noisy and gate on more homogeneous populations. Also, I do an additional doublet discrimination step of SSCH vs SSCW in my manual analysis, do you think this would be helpful here?
