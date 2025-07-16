lib: circle.f90 clean
	ifx -shared -fpic -O0 -g circle.f90 -o libclass_circle.so
	nm ./libclass_circle.so > nm_libclass_circle.so.txt

exe: circle_test.f90 lib
	ifx -L. -lclass_circle -o exe circle_test.f90 

clean:
	rm -f ./exe
	rm -f ./*.so

run: exe lib
	# ./exe
	python t.py


dbg: lib
	python -m pdb t.py