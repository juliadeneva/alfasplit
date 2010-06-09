alfasplit: alfasplit.c 
	gcc -g -o alfasplit alfasplit.c alfaoff.c azza.c prec.c radec.c deg_trig.c cal2mjd.c y.tab.c head_parse.c -I/home/cole/deneva/alfasplit -lm 
	cp -f ./alfasplit ~/bin64/.

alfabeam: alfabeam.c
	gcc -g -o alfabeam alfabeam.c alfaoff.c azza.c prec.c radec.c deg_trig.c cal2mjd.c -I/home/cole/deneva/alfasplit -lm 

clean:
	rm alfasplit alfabeam *.o



