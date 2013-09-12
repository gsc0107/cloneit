.PHONY:all clean



all:allenz.js pbsrf2.js pgbt9.js

allenz.js:
	curl "ftp://ftp.neb.com/pub/rebase/allenz.txt" |\
	awk '/^<1>/ { N=substr($$0,4);next;} /^<5>/{S=substr($$0,4);next;} /^<7>/{C=length(substr($$0,4)); if(C==0 || index(S,"^")==0) next; printf("%s\t%s\t%d\n",N,S,C);next;}' |\
	LC_ALL=C sort -t '	' -k2,2 -k3,3n -u |\
	LC_ALL=C sort -t '	' -k1,1 |\
	awk -F '	' 'BEGIN {printf("function makeRebase()\n{\nvar r=new Rebase();\n");} {printf("r.add(\"%s\",\"%s\",%s);\n",$$1,$$2,$$3);} END { printf("return r;\n}\n");}' > $@

pbsrf2.js: ../test/fasta/pbs_RF2.fa
	echo -n "var pbsrf2=Sequence.create('pBSRF2','" > $@ && \
	grep -v '>' $< | tr -d '\n' >> $@ && \
	echo  "');" >> $@


pgbt9.js: ../test/fasta/pGBT9.fa
	echo -n "var pgbt9=Sequence.create('pGBT9','" > $@ && \
	grep -v '>' $< | tr -d '\n' >> $@ && \
	echo  "');" >> $@


clean:
	rm -f allenz.js
