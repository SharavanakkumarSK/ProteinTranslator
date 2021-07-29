# ProteinTranslator 

def transcription(dna_seq):
    no_of_dna_bases = len(dna_seq)
    print("No of DNA Bases:", no_of_dna_bases)
    for i in dna_seq:
        rna_seq = dna_seq.replace('T', 'U')
        return "RNA Sequence: " +  rna_seq
    
def translation(dna_seq):
    no_of_dna_bases = len(dna_seq)
    print("No of RNA Bases:", no_of_dna_bases)
       
    genetic_code = {"GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
                    "TGC":"C", "TGT":"C", "GAC":"D", "GAT":"D",
                    "GAA":"E", "GAG":"E", "TTC":"F", "TTT":"F",
                    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
                    "CAC":"H", "CAT":"H", "ATA":"I", "ATC":"I",
                    "ATT":"I", "AAA":"K", "AAG":"K", "TTA":"L",
                    "TTG":"L", "CTA":"L", "CTC":"L", "CTG":"L", 
                    "CTT":"L", "ATG":"M", "AAC":"N", "AAT":"N",
                    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
                    "CAA":"Q", "CAG":"Q", "AGA":"R", "AGG":"R", 
                    "CGA":"R", "CGC":"R", "CGT":"R", "CGG":"R", 
                    "AGC":"S", "AGT":"S", "TCA":"S", "TCC":"S", 
                    "TCG":"S", "TCT":"S", "ACA":"T", "ACC":"T", 
                    "ACG":"T", "ACT":"T", "GTA":"V", "GTC":"V",
                    "GTG":"V", "GTT":"V", "TGG":"W", "TAC":"Y", 
                    "TAT":"Y", "TAG":"!", "TAA":"!!", "TGA":"!!!"}
    protein_seq = ""
    for b in range(no_of_dna_bases):
        if(dna_seq[b : b + 3] == "ATG"):
            print("Start Codon Founded!")
            print("Position of Start Codon -", b)
                   
            for j in range(b, no_of_dna_bases, 3):
                amino_acid = genetic_code[dna_seq[j : j + 3]]
                if(amino_acid == "!"):
                    print("Stop Codon Founded! - UAG")
                    print("Position of Stop Codon -", j)
                    break
                elif(amino_acid == "!!"):
                    print("Stop Codon Founded! - UAA")
                    print("Position of Stop Codon -", j)
                    break
                elif(amino_acid == "!!!"):
                    print("Stop Codon Founded! - UGA")
                    print("Position of Stop Codon -", j)
                    break
                else:
                    protein_seq = protein_seq + amino_acid
            print("Amino Aicd Sequence: ", protein_seq)
            no_amino_acids = len(protein_seq)
            print("No of Amino Acids:", no_amino_acids)
            
# Type - 1 input, for entering the sequence in a query box:

dna_seq = str(input("DNA Sequence: ")).upper()
result_tscn = transcription(dna_seq)
print(result_tscn)
result_tsln = translation(dna_seq)
print(result_tsln)

# Type - 2 input, for opening the sequence in a text file:

inp_seq = open("'_directoryname_':\'_filename_'.txt", "r")
dna_seq = inp_seq.read()
result_tscn = transcription(dna_seq)
print(result_tscn)
result_tsln = translation(dna_seq)
print(result_tsln)
