### a test-case for the ELPA development. 

The original files are in `pasta-test.tar.gz`

Idn bijlage een voorbeeldje: 3x3x3 spherische clusters op relatief hoge
dichtheid. Er zijn vier bestanden

  * run*sh: het runscriptje; geconfigureerd voor de setup van Nikolai,
    maar het zou je duidelijk moeten maken hoe de andere bestandjes aan
    MOCCa gevoederd worden.
  * data* : de input data voor MOCCa
  * out*  : dit is de STDOUT voor de code; die bevat eigenlijk alle
    fysieke informatie.
  * inp*pot: dit bevat onze initialisatie voor de code, bepaald door een
    semiklassieke benadering van het probleem.

Nog wat opmerkingen:

  * Zoals je kan zien in het runscript, kan dit voorbeeld runnen op 8
    nodes; Nikolai heeft het ook getest op 6 & 10 nodes.
  * De rekentijd is wel serieus als je evenveel iteraties wilt doen: ~9u
    voor 400 iteraties. Je kan echter zoveel iteraties doen als je
    budget hebt door met het keyword "maxiter" in het data* bestandje te
    spelen. Schatting: ~1.5 minuten per iteratie op 6 nodes.

The test needs to be run with Tantalus.BXL.exe.
