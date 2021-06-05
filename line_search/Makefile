#Checar en http://www.cs.dartmouth.edu/~campbell/cs50/buildlib.html
#para generacion de una biblioteca estatica a partir de los
#archivos .o generados a partir de la compilacion de un archivos .cpp

#Checar en https://www.gnu.org/software/make/manual/html_node/Wildcard-Function.html
#para analizar el funcionamiento de las funciones "patsubst" y "wildcard"

CCOMP=gcc		#compilador
CLINK=gcc		#linkeo
LFLAGS=-std=c99 	#bandera de debuggeo
CFLAGS=-std=c99	-c	#banderas de compilacion
EXEC=main.exe
PDFDIR=docs
DOCDIR:=latex
OBJDIR:=obj
CFILEDIR:=src
HFILEDIR:=src
BINDIR:=bin
DOCS:=$(wildcard $(PDFDIR)/*.*)
README:=$(wildcard $(DOCDIR)/*.tex)
CFILES:=$(wildcard $(CFILEDIR)/*.cpp)
LIBS:=$(wildcard $(HFILEDIR)/*.hpp)
OBJS:=$(patsubst $(CFILEDIR)/%.cpp,$(OBJDIR)/%.o,$(CFILES))
PDFS:=$(patsubst $(DOCDIR)/%.tex,$(PDFDIR)/%.pdf,$(README))
DIREXEC:=$(addprefix $(BINDIR)/,$(EXEC))

$(DIREXEC): $(OBJS) $(LIBS)
	$(CCOMP) $(LFLAGS) -o $@ $(OBJS)

$(OBJDIR)/%.o: $(CFILEDIR)/%.cpp $(LIBS) 
	$(CCOMP) $(CFLAGS) $< -o $@

optimized:
	make CFLAGS="-O2 -Wall -std=c99 -c"; LFLAGS="-O2 -Wall"

debugging:
	make CFLAGS="-g -std=c99 -c"; LFLAGS="-g -std=c99 -o"

clean:
	rm $(OBJS) $(DIREXEC) $(DOCS)

readme: $(PDFS)
	@echo $(PDFS)

$(PDFDIR)\%.pdf: $(DOCDIR)\%.tex
	pdflatex -output-directory=.\docs $<

run:
	./$(DIREXEC) ${args}

test:
	make run args="Cholesky data/A.dat data/b.dat"
	
help:
	@echo 	"Archivos a compilar: obj/main.o  obj/listas.o  obj/memo.o y bin/ejecutable\n Opciones de compilacion: make, make optimized, make debugging"