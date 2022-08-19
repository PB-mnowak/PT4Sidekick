protein_db = {
    'ss': [
        'ATGTACAGGATGCAACTCCTGTCTTGCATTGCACTAAGTCTTGCACTTGTCACGAATTCA',
        'ATGTACAGGATGCAACTCCTGTCTTGCATTGCACTAAGTCTTGCACTTGTCACGAATTCG',
        
    ],
    'end': [
        'ACAAAGAGCTTCAACAGGGGAGAGTGTTAG', # Light chain - 
        'AAGAGCCTCTCCCTGTCTCCGGGTAAATGA', # Heavy chain - IgG1
        'GCCTCTCCCTGTCTCCGGGTAAATAAAC'    # Heavy chain - IgG4
    ]
}

lg_units = {
    'Kg': 1,
    'g': 2,
    'mg': 3,
    'μg': 4,
    'ng': 5,
    'L': 6,
    'mL': 7,
    'μL': 8,
    'mg/mL': 9,
    'g/L': 10,
    'M': 11,
    'mM': 12,
    'μM': 13,
    'μg/mL': 14,
    'ng/μL': 15,
    'ratio': 16,
    '%': 17,
    'g/mL': 18,
}

lg_units_conversion = {
    6: 10**6,
    7: 10**3,
    8: 1,
}

conc_units = {
    'kg': 10**12,
    'g': 10**9,
    'mg': 10**6,
    'μg': 10**3,
    'ng': 1,
    'pg': 10**-3,
    'l': 10**6,
    'ml': 10**3,
    'μl': 1,
    'nl': 10**-3,
    'pl': 10**-6,
}

box_types = {
    "Primer": 0,
    "General": -1,
    "Tube": -1,
    "Seed": -2,
    "Fly": -3,
    "Bacterium": -4,
    "Bacteria": -4,
    "CellLine": -5,
    "Cell Line": -5,
    "Cellline": -5,
    "Tissue": -6,
    "Antibody": -7,
    "Plasmid": -8,
    "Glycerol": -9,
    "Enzyme": -10,
    "Consumable": -11,
    "Yeast": -12,
    "Fungus": -13,
    "Virus": -14,
    "Protein": -15,
    "Lipid": -16,
    "Worm": -17,
    "Sequence": -18,
    "Zebrafish": -19,
    "Gene": -20,
    "Compound": -21,
}

sheet_verification = {
    'sheets': {
        'Plasmids',
        'Proteins'
        },
    'Plasmids': {'ID',
                 'Stock name (Plasmid)',
                 'Filtered',
                 'Conc',
                 'Volume',
                 'Box',
                 'Position',
                 'Link - Stock',
                 'SysID',
                 'Plasmid inventory name',
                 'Link - Inv',
                 'Transfer',
                 'New Box',
                 'New Position'
                 },
    'Proteins': {
        'Plasmid SysID',
        'Plasmid inventory name',
        'POI SysID',
        'POI ID',
        'POI name',
        'Concentration',
        'Tag',
        'MW',
        'Length',
        'pI',
        'Purification method',
        'A0.1% (Ox)',
        'A0.1% (Red)',
        'Description',
        'Stock 1',
        'Stock 2',
        'Stock 3',
        'Stock 4',
        'Stock 5'
        },  
}