XLF= ifort -O3 -static -heap-arrays 1024

bhf.exe :   bhf-modules.o bhf-main.o bhf-library.o bhf-matrix.o  bhf-closedshells.o heff-diagrams.o  
	       ${XLF} -o  bhf.exe  bhf-modules.o bhf-main.o bhf-library.o bhf-matrix.o bhf-closedshells.o heff-diagrams.o

bhf-modules.o :   bhf-modules.f90
		${XLF} -c  bhf-modules.f90


bhf-main.o :   bhf-main.f90
		${XLF} -c bhf-main.f90

bhf-closedshells.o :   bhf-closedshells.f90
		${XLF} -c  bhf-closedshells.f90


heff-diagrams.o :   heff-diagrams.f90
		${XLF} -c  heff-diagrams.f90

bhf-library.o :   bhf-library.f
	          ${XLF} -c bhf-library.f 

bhf-matrix.o :   bhf-matrix.f90
		${XLF} -c   bhf-matrix.f90
