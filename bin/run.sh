gnuplot draw.gp

# latex
latex  Case-xx-SOD-rho.tex
latex  Case-xx-SOD-u.tex
latex  Case-xx-SOD-p.tex
latex  Case-xx-SOD-e.tex

# ps
dvips Case-xx-SOD-rho.dvi
dvips Case-xx-SOD-u.dvi
dvips Case-xx-SOD-p.dvi
dvips Case-xx-SOD-e.dvi

# eps
ps2eps -f Case-xx-SOD-rho.ps
ps2eps -f Case-xx-SOD-u.ps
ps2eps -f Case-xx-SOD-p.ps
ps2eps -f Case-xx-SOD-e.ps

# dvipdf 
#dvipdf Case-xx-SOD-rho.dvi
#dvipdf Case-xx-SOD-u.dvi
#dvipdf Case-xx-SOD-p.dvi
#dvipdf Case-xx-SOD-e.dvi

rm *.tex *.aux *.dvi *.log *-inc.eps *.ps 
