frame_dict = {"TTT":"F|Phe", "TTC":"F|Phe", "TTA":"L|Leu", "TTG":"L|Leu", "TCT":"S|Ser", "TCC":"S|Ser", "TCA":"S|Ser",
              "TCG":"S|Ser", "TAT":"Y|Tyr", "TAC":"Y|Tyr", "TAA":"*|Stp", "TAG":"*|Stp", "TGT":"C|Cys", "TGC":"C|Cys",
              "TGA":"*|Stp", "TGG":"W|Trp", "CTT":"L|Leu", "CTC":"L|Leu", "CTA":"L|Leu", "CTG":"L|Leu", "CCT":"P|Pro",
              "CCC":"P|Pro", "CCA":"P|Pro", "CCG":"P|Pro", "CAT":"H|His", "CAC":"H|His", "CAA":"Q|Gln", "CAG":"Q|Gln",
              "CGT":"R|Arg", "CGC":"R|Arg", "CGA":"R|Arg", "CGG":"R|Arg", "ATT":"I|Ile", "ATC":"I|Ile", "ATA":"I|Ile",
              "ATG":"M|Met", "ACT":"T|Thr", "ACC":"T|Thr", "ACA":"T|Thr", "ACG":"T|Thr", "AAT":"N|Asn", "AAC":"N|Asn",
              "AAA":"K|Lys", "AAG":"K|Lys", "AGT":"S|Ser", "AGC":"S|Ser", "AGA":"R|Arg", "AGG":"R|Arg", "GTT":"V|Val",
              "GTC":"V|Val", "GTA":"V|Val", "GTG":"V|Val", "GCT":"A|Ala", "GCC":"A|Ala", "GCA":"A|Ala", "GCG":"A|Ala",
              "GAT":"D|Asp", "GAC":"D|Asp", "GAA":"E|Glu", "GAG":"E|Glu", "GGT":"G|Gly", "GGC":"G|Gly", "GGA":"G|Gly",
              "GGG":"G|Gly"
              }
# writes the file
def writefile(col):
    DNA, revDNA, frames = col
    with open("6frame.txt", "w") as file:
        file.write(f"3+   {''.join(frame_dict[codon].split('|')[1] for codon in frames[3] if codon in frame_dict.keys())}\n"
                   f"2+  {''.join(frame_dict[codon].split('|')[1] for codon in frames[2] if codon in frame_dict.keys())}\n"
                   f"1+ {''.join(frame_dict[codon].split('|')[1] for codon in frames[1] if codon in frame_dict.keys())}\n"
                   f"   {DNA}\n"
                   f"   {revDNA}\n"
                   f"1- {''.join(frame_dict[codon].split('|')[1] for codon in frames[-1] if codon in frame_dict.keys())}\n"
                   f"2-  {''.join(frame_dict[codon].split('|')[1] for codon in frames[-2] if codon in frame_dict.keys())}\n"
                   f"3-   {''.join(frame_dict[codon].split('|')[1] for codon in frames[-3] if codon in frame_dict.keys())}\n\n"
                   
                   f"3+   {' '.join(frame_dict[codon].split('|')[0] for codon in frames[3] if codon in frame_dict.keys())}\n"
                   f"2+  {' '.join(frame_dict[codon].split('|')[0] for codon in frames[2] if codon in frame_dict.keys())}\n"
                   f"1+ {' '.join(frame_dict[codon].split('|')[0] for codon in frames[1] if codon in frame_dict.keys())}\n"
                   f"   {DNA}\n"
                   f"   {revDNA}\n"
                   f"1- {' '.join(frame_dict[codon].split('|')[0] for codon in frames[-1] if codon in frame_dict.keys())}\n"
                   f"2-  {' '.join(frame_dict[codon].split('|')[0] for codon in frames[-2] if codon in frame_dict.keys())}\n"
                   f"3-   {' '.join(frame_dict[codon].split('|')[0] for codon in frames[-3] if codon in frame_dict.keys())}")
# returns reverse complement
def revcomp(strand):
    newdna = ''

    for letter in strand:
        if letter == "A":
            newdna += "T"
        elif letter == "C":
            newdna += "G"
        elif letter == "G":
            newdna += "C"
        else:
            newdna += "A"

    return newdna

# generates a 6frame translation
def sixframe(strand):
    revstrand = revcomp(strand[:-1])
    frame = {
        1: [],
        2: [],
        3: [],
        -1: [],
        -2: [],
        -3: [],
    }
    for i in range(len(strand)):
        frame[(i%3)+1].append(strand[i:i+3])
        frame[(i%-3)-1].append(revstrand[i:i+3])

    return strand, revstrand, frame



# Reads DNA.txt and uses 6frame() to generate
def fileread():
    with open("dna.txt", "r") as file:
        return  file.read().strip().upper()


if __name__ == ("__main__"):
    k = 0

    while k != any([1, 2]):
        k = int(input("wil je handmatig een DNA strand invoeren(1)\n"
                      "of het bestand DNA.txt inlezen(2)\n"))
        if k == 1:
            DNAlegal = False
            DNA = input("voer een DNA string in:\n").upper()
            if all(char.isalpha() and char in "ATGC" for char in DNA):
                DNAlegal = True

            while not DNAlegal:
                DNA = input("foute tekens gevonden, voer een nieuwe DNA string in:\n")
                if all(char.isalpha() and char in "ATGC" for char in DNA):
                    DNAlegal = True

            writefile(sixframe(DNA))
        else:
            writefile(sixframe(fileread()))
