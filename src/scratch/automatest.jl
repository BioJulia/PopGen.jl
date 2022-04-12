# automa for genepop parsing

import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp;

machine = let
    # Define genepop syntax
    header = re"[^\n]*\n"
    popsep = re"[PpOo]+" * re"[\r\n]"
    sampleid = re"[a-zA-Z0-9_\-.]+"
    samplesep = re",[ \t]*"
    loci = re.rep(re"[a-zA-Z0-9_\-.]+[\r\n,]*")
    genotypes = re.rep(re"[\-]?[0-9]+[ |\t\r\n]")
    record = sampleid * samplesep * genotypes 
    genepop = header * loci * popsep * re.rep(record)

    # Compile the regex to a FSM
    Automa.compile(genepop)
end;
