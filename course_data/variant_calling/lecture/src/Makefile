# notes=1 make
#
# ln -rs ~/wtxt/slides/pandoc/pandoc.Makefile Makefile
# ln -rs ~/wtxt/slides/pandoc/pandoc.style.tex
# ln -rs ~/wtxt/slides/pandoc/pandoc.style.notes.tex
# ln -rs ~/wtxt/slides/pandoc/pandoc.filter
# ln -rs ~/wtxt/slides/pandoc/pandoc.includes
#

all:
	@./pandoc.includes -m > rmme.mak && make -f rmme.mak && rm -f rmme.mak

tex:
	@./pandoc.includes -m > rmme.mak && make -f rmme.mak tex && rm -f rmme.mak

clean:
	@rm -f *.log *.aux *.out *.nav *.snm *.toc && ./pandoc.includes -m > rmme.mak && make -f rmme.mak clean && rm -f rmme.mak

