lib: elpa_object_interface.f90 
	ifx -shared -fpic \
		-I/apps/antwerpen/skylake/rocky9/ELPA/2024.05.001-intel-2024a/include/elpa-2024.05.001/modules/ \
		elpa_object_interface.f90 -o libelpa_object_interface.so
	nm ./libelpa_object_interface.so > nm_libelpa_object_interface.so.txt

clean:
	rm ./*.so
