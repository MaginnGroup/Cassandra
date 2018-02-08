latex user_guide.tex
latex user_guide.tex
bibtex user_guide
bibtex user_guide
latex user_guide.tex
dvips -o user_guide.ps user_guide.dvi

#ps2pdf user_guide.ps

# For MacOS
dvipdfmx user_guide.dvi
