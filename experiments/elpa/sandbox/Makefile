lib: circle.f90 clean
	ifx -shared -fpic circle.f90 -o libclass_circle.so
	nm ./libclass_circle.so > nm_libclass_circle.so.txt

exe: circle_test.f90 libclass_circle.so
	ifx -L. -lclass_circle -o exe circle_test.f90 

clean:
	rm -f ./exe
	rm -f ./*.so

run: exe lib
	# ./exe
	python pycircle_test.py