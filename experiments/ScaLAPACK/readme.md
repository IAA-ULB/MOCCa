following [this link](https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack) to get started with scalapack.

 # [ScaLAPACK tutorial on Netlib](https://www.netlib.org/scalapack/tutorial)

 ## [slide 22](https://www.netlib.org/scalapack/tutorial/sld022.htm)

 each global data object is assigned an array descriptor, which 
 - contains information required to establish mapping between a global array entry and its corresponding process and memory location
 - id differentiated by de DTYPE_(first entry) in the descriptor
 - provides a flexible framework to eqsily specify additional data distributions or matrix types

[ScaLAPACK Tutorial (pdf) by Jack Dongarra and Susan Blackford](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=7b2ae178767a55906c7145221b033f1c49ec5ab3)

block cyclic data distribution: 

    verdeel de globale matrix in blokken 
    zet de grid-rij index op 1
    loop over blok-rijen
        zet de grid-kolom index op 1
        loop over blok-kolommen in huidige blok-rij
            geef het blok aan het huidige grid-proces
            increment grid-kolom index
        increment grid-row index

![een 5x5 matrix gepartitioneerd in 2x2 blokken verdeeld over een 2x2 grid](image.png)

[PRACE ScaLAPACK Tutorial](https://events.prace-ri.eu/event/1286/attachments/1667/3912/ScaLAPACK_PTC.pdf)

[Python ScaLAPACK wrapper](https://pypi.org/project/PyScalapack/)
