import requests
import xml.etree.ElementTree as ET
import os

class Protein():
    
    def __init__(self) -> None:
        self.aa_dict = {'A': ('Ala', 'Alanine', 71.0788), 
                        'R': ('Arg', 'Arginine', 156.1875), 
                        'N': ('Asn', 'Asparagine', 114.1038), 
                        'D': ('Asp', 'Aspartic Acid', 115.0886), 
                        'C': ('Cys', 'Cysteine', 103.1388), 
                        'Q': ('Gln', 'Glutamine', 128.1307), 
                        'E': ('Glu', 'Glutamic Acid', 129.1155), 
                        'G': ('Gly', 'Glycine', 57.0519), 
                        'H': ('His', 'Histidine', 137.1411), 
                        'I': ('Ile', 'Isoleucine', 113.1594), 
                        'L': ('Leu', 'Leucine', 113.1594), 
                        'K': ('Lys', 'Lysine', 128.1741), 
                        'M': ('Met', 'Methionine', 131.1926), 
                        'F': ('Phe', 'Phenylalanine', 147.1766), 
                        'P': ('Pro', 'Proline', 97.1167), 
                        'S': ('Ser', 'Serine', 87.0782), 
                        'T': ('Thr', 'Threonine', 101.1051), 
                        'W': ('Trp', 'Tryptophan', 186.2132), 
                        'Y': ('Tyr', 'Tyrosine', 163.176), 
                        'V': ('Val', 'Valine', 99.1326)}

        self.glyco_sites = []
        self.ds_sites = []
        self.lipid_sites = []
        
        pass

    def uniprot_data(self, up_id=None):   
        
        if up_id == None:
            up_id = input("Enter UniProt ID: ")
        
        link_backbone = "https://www.uniprot.org/uniprot/"
        link_format = ".xml"
        url_str = link_backbone + up_id + link_format

        up_xml = requests.get(url_str)
        tree = ET.fromstring(up_xml.text)
        

        
        self.up_seq = tree.find('{http://uniprot.org/uniprot}sequence')
        
        for feature in tree.iter('{http://uniprot.org/uniprot}feature'):
            if "glycosylation site" in feature.attrib["type"]:
                self.glyco_sites.append(int(feature[0][0].attrib["position"]))
            elif "disulfide bond" in feature.attrib["type"]:
                begin = int((feature[0][0].attrib["position"]))
                end = int((feature[0][1].attrib["position"]))
                self.ds_sites.append((begin, end))
            elif "lipid moiety-binding region" in feature.attrib["type"]:
                self.lipid_sites.append(int(feature[0][0].attrib["position"]))
        
        # print(self.up_seq)
        # print(self.glyco_sites)
        # print(self.ds_sites)
        # print(self.lipid_sites)
        # page = requests.get(link_backbone + up_id + link_format)

    # Generate list of aa (list) and calculate it's length (int)
    def add_sequence(self, seq=None):
        if not seq:
            seq = input("Enter protein sequence: ").strip().upper()
        self.seq = "".join([aa for aa in seq if aa in self.aa_dict])
        self.seq_len = len(self.seq)

    def seq_comp(self):
        pass

    # Calculate aa distribution (dict)
    def aa_distr(self):
        self.aa_count = dict((aa, self.seq.count(aa)) for aa in self.aa_dict)
        self.aa_proc = dict((aa, self.aa_count[aa] / self.seq_len) for aa in self.aa_dict)

    def dipeptides(self):
        
        dp_dict = {}
        i = 0
        
        for aa in self.seq:
            dp = self.seq[i:i+2]
            dp_num = dp_dict.setdefault(dp, 0)
            dp_dict[dp] = dp_num + 1
            i += 1
            
        self.dp_count = dp_dict

    # Calculate protein mass
    def prot_mass(self):
        self.mass = sum(self.aa_dict.get(aa, 0)[2] for aa in self.seq) + 18.01528

    # Estimate protein absorption coefficients
    def abs_coeff(self):
        base = 5500 * self.aa_count["W"] + 1490 * self.aa_count["Y"]
        
        self.e_red = base
        self.e_ox = base + 125 * (self.aa_count["C"] // 2)
        self.a_red = base / self.mass
        self.a_ox = self.e_ox / self.mass

    def pI(self):
        
        # Calculates charge of aa at given pH based on pKa
        def charge(ph, pka):
            return 1 / (10 ** (ph - pka) + 1)
        
        # def net_charge():
        #     return
        
        # Set range of pI search and start pH
        pH_min, pH, pH_max = 2, 8, 14.0
        
        base_pka = {'K': 10.0,
                    'R': 12.0,
                    'H': 5.98}
        acid_pka = {'D': 4.05,
                    'E': 4.45,
                    'C': 9.0,
                    'Y': 10.0}
        n_term_dict = {'A': 7.59, 
                       'M': 7.0, 
                       'S': 6.93, 
                       'P': 8.36, 
                       'T': 6.82, 
                       'V': 7.44, 
                       'E': 7.7}
        c_term_dict = {'D': 4.55, 
                       'E': 4.75}

        n_pka = n_term_dict.get(self.seq[0], 7.7)
        c_pka = c_term_dict.get(self.seq[-1], 3.55)
        
        p_charge = sum([charge(pH, base_pka[aa]) * self.aa_count[aa] for aa in base_pka.keys()]) + charge(pH, n_pka)
        n_charge = sum([charge(acid_pka[aa], pH) * self.aa_count[aa] for aa in acid_pka.keys()]) + charge(c_pka, pH)
        net_charge = p_charge - n_charge

        while abs(net_charge) > 0.0001:
            if  net_charge > 0:
                pH_min = pH
            else:
                pH_max = pH
            
            pH = (pH_min + pH_max) / 2
            
            p_charge = sum([charge(pH, base_pka[aa]) * self.aa_count[aa] for aa in base_pka.keys()]) + charge(pH, n_pka)
            n_charge = sum([charge(acid_pka[aa], pH) * self.aa_count[aa] for aa in acid_pka.keys()]) + charge(c_pka, pH)
            net_charge = p_charge - n_charge

        self.pI = pH

    # Print calculated results
    def show_analysis(self):
        
        print()
        
        self.row_len = 50
        self.delim = 10
        
        print(f"Sequence:", end="")
        seq_cut = [(i, self.seq[i:i+self.delim]) for i in range(0, self.seq_len, self.delim)]
        for n in seq_cut:
            if n[0] % self.row_len != 0:
                print(n[1], end=" ")
            else:
                print(f"\n{n[0]:>4d}\t{n[1]}", end=" ")
        print("\n")

        print(f"Amino acids:")
        for aa in self.aa_count:
            print(f"{aa} | {self.aa_dict[aa][0]} \t{self.aa_count[aa]:>2d}\t{self.aa_proc[aa]:>4.1f}%")
            
        print()
        
        print(f"Length:\t\t{self.seq_len}")
        
        print(f"pI:\t\t{self.pI:.2f}")
        
        print(f"Mass:\t\t{self.mass:.2f} Da")
        
        print(f"E(ox):\t\t{self.e_ox:d}")
        
        print(f"E(red):\t\t{self.e_red:d}")
        
        print(f"A0.1%(ox):\t{self.a_ox:.3f}")
        
        print(f"A0.1%(red):\t{self.a_red:.3f}")

        print()
        if self.ds_sites:
            print(f"Disulfide bonds: [{len(self.ds_sites)}]" + "\t" + "    v - v    ")
            for ds1, ds2 in self.ds_sites:
                print(f"\t{ds1:<3d} - {ds2:<3d}\t{self.seq[ds1-5:ds1]} - {self.seq[ds2-1:ds2+4]}")
            print()
                
        if self.glyco_sites:
            print(f"Glycosylation: [{len(self.glyco_sites)}]" + "\t" + "  v  ")
            for g in self.glyco_sites:
                print(f"\t{g}\t{self.seq[g-1]}\t{self.seq[g-3:g+2]}")
            print()
                
        if self.lipid_sites:
            print(f"Lipidation: [{len(self.lipid_sites)}]" +"\t" * 2 + "  v  ")
            for l in self.lipid_sites:
                print(f"\t{l}\t{self.seq[l-1]}\t{self.seq[l-3:l+2]}")
            print()


def run_script():
    os.system('python3 "P:\_research group folders\PT Proteins\protein.py"')

def load_file(address):
    with open(address, "r") as f:
        for line in f.readline():
            if ">" in line:
                up_id = line.split()

# if __name__ == "__main__":
#     # up_id = "P41221"
#     # seq = "MKKSIGILSPGVALGMAGSAMSSKFFLVALAIFFSFAQVVIEANSWWSLGMNNPVQMSEVYIIGAQPLCSQLAGLSQGQKKLCHLYQDHMQYIGEGAKTGIKECQYQFRHRRWNCSTVDNTSVFGRVMQIGSRETAFTYAVSAAGVVNAMSRACREGELSTCGCSRAARPKDLPRDWLWGGCGDNIDYGYRFAKEFVDARERERIHAKGSYESARILMNLHNNEAGRRTVYNLADVACKCHGVSGSCSLKTCWLQLADFRKVGDALKEKYDSAAAMRLNSRGKLVQVNSRFNSPTTQDLVYIDPSPDYCVRNESTGSLGTQGRLCNKTSEGMDGCELMCCGRGYDQFKTVQTERCHCKFHWCCYVKCKKCTEIVDQFVCK"
    
#     poi = Protein()
#     # poi.uniprot_data(up_id=up_id)
#     # poi.add_sequence(seq=seq)
#     poi.add_sequence()
#     poi.aa_distr()
#     poi.prot_mass()
#     poi.abs_coeff()
#     poi.pI()
#     poi.show_analysis()
