######################################################
#
# makefile for cloneit 
# cgi version 2.1
#
######################################################

###############################
#
# path on DISK
#
#   to binary executable
#
BIN_PATH=/home/apache2/cgi-bin/caltech/
#
#   to directory containing html files, allenz etc...
#
DATAS_PATH=\"/usr/local/src/cloneitcgi/\" 
#
#   to html files
#
HELP_PATH_DISK=/home/apache2/htdocs/saf_manuals/cloneit/
###############################
#
# URLS 
#
#  to cgi bin program
#
CGI_PATH=/cgi-bin/caltech/
#
#   to html file files
#
HELP_PATH=/saf_manuals/cloneit/
###########################################################
CC=gcc
CFLAGS= -O3 -Wall
BIGSED=	sed 's/__MY_MAIL__/lindenb@cng.fr/g' | \
		sed 's%__BIBLIO__%http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve\&db=PubMed\&list_uids=9682060\&dopt=Abstract%g' |\
		sed 's%__BMR__%http://locus.jouy.inra.fr/virim/L002.htm%g' |\
		sed 's%__REBASE__%http://rebase.neb.com%g' |\
		sed 's%__CGI_PATH__%${CGI_PATH}%g' |\
		sed "s%__HELP_PATH__%${HELP_PATH}%" 


all: cloneitcgi classic_allenz html_files

install:
	cp cloneitcgi ${BIN_PATH}/cloneitcgi
	cp *.html     ${HELP_PATH_DISK}

cloneitcgi:CloneItCGI.c html_files classic_allenz
	${CC} ${CFLAGS} -DDATAS_PATH="${DATAS_PATH}" -o $@ CloneItCGI.c

archive: 
	make cloneitcgi.tar.gz

cloneitcgi.tar.gz:
	-rm -f cloneitcgi.tar.gz 
	tar -cvf cloneitcgi.tar CloneItCGI.c Error.html_ Polylinkers.Set README allenz  cloneit.jpeg cloneit_cgi0.html_ cloneit_cgi1.html_ cloneit_cgi2.html_ \
		cloneithelp.html_ head.html_ makefile misc.html_ pBS_RF2.txt pGAD424.txt \
	    tail.html_
	gzip cloneitcgi.tar

classic_allenz: allenz
	egrep "^<[0-8]>" allenz | \
	sed 's/^<1>/#<1>/' | sed 's/^</=</' | tr -d "\n" | tr -s "#" "\n" | \
	egrep "(<1>EcoRI=|<1>BamHI=|<1>HindIII=|<1>PstI=|<1>SmaI=|<1>XmaI=|<1>AvrII=|<1>BglII=|<1>KpnI=|<1>NheI=|<1>NcoI=|<1>NotI=|<1>PvuII=|<1>SalI=|<1>XhoI=|<1>XbaI=)" | \
	tr -s "=" "\n" > $@

html_files: head.html tail.html misc.html CloneIt_Help.html cloneit_cgi.html example.html Error.html

Error.html: Error.html_
	cat  Error.html_  | ${BIGSED}  > $@    

head.html: head.html_
	cat head.html_ | ${BIGSED} | sed "s/__DATE__/`date`/" > $@


tail.html: tail.html_
	cat tail.html_ | ${BIGSED} > $@

misc.html: misc.html_
	cat misc.html_  | ${BIGSED} > $@

CloneIt_Help.html: head.html cloneithelp.html_ misc.html tail.html
	cat head.html cloneithelp.html_ misc.html tail.html  |\
	sed 's%Online</TITLE>%Help</TITLE>%g' | ${BIGSED} > $@


cloneit_cgi.html: cloneit_cgi0.html_  cloneit_cgi1.html_ cloneit_cgi2.html_ misc.html
	cat  cloneit_cgi0.html_ cloneit_cgi1.html_ cloneit_cgi2.html_ misc.html  |\
	${BIGSED} > $@

example.html:  head.html cloneit_cgi0.html_  cloneit_cgi1.html_ cloneit_cgi2.html_ misc.html  tail.html pBS_RF2.txt pGAD424.txt
	cat head.html cloneit_cgi0.html_ pGAD424.txt cloneit_cgi1.html_ pBS_RF2.txt cloneit_cgi2.html_  misc.html tail.html  |\
	sed 's/Value="vector" name="V_name"/Value="pGAD424" name="V_name"/' |\
	sed 's/Value="insert" name="I_name"/Value="PBS-RF2" name="I_name"/' |\
	sed 's%^<form name%<p><b>Let s clone the NSP5(RF2) Rotavirus protein <b>in frame</b> into the pGAD424 vector (From CLONTECH)</b></p><form name%' |\
	sed '/frame5/s/value="1"/value="1" chekcked/' |\
	sed 's%Online</TITLE>%Example</TITLE>%g' | ${BIGSED} > $@

clean:
	rm -f *.o ${BIN_PATH}/cloneitcgi cloneitcgi head.html classic_allenz tail.html \
		misc.html CloneIt_Help.html cloneit_cgi.html example.html Error.html
