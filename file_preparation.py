import os
import csv
from typing import Set, List, Dict
from collections import defaultdict


def extract_Entrez(input_path: str, output_path: str, taxon_ids: List[str]):
    """
    Extract species specific Entrez identifiers.
    :param input_path: path to the input file with Entrez identifiers
    :param output_path: path to the output file for species specific Entrez identifiers
    :param taxon_ids: list of numerical taxon IDs (NCBI) for which to extract the identifiers
    """
    with open(input_path, 'r', encoding='utf-8') as in_file, open(output_path, 'w', encoding='utf-8') as out_file:
        out_file.write('Entrez ID\tGene name\tEntrez Synonyms\n')
        for line in csv.DictReader(in_file, delimiter='\t'):
            # skip entries associated with other species
            if line['#tax_id'] not in taxon_ids:
                continue
            # skip lines without Entrez ID
            if not line['GeneID'] or line['GeneID'] == '-':
                continue
            name = '' if line['Symbol'] == '-' else line['Symbol']                                  # type: str
            synonyms = '' if line['Synonyms'] == '-' else line['Synonyms'].replace('|', ', ')       # type: str
            out_file.write('\t'.join([line['GeneID'], name, synonyms]) + '\n')


def extract_ORegAnno(input_path: str, output_path: str, species: List[str]):
    """
    Extract species specific regulatory information from an ORegAnno data file.
    :param input_path: path to the input file with all annotations
    :param output_path: path to the output for species specific regulatory annotations
    :param species: list of species names for which to extract regulatory annotations
    """
    with open(input_path, 'r') as original_file, open(output_path, 'w') as new_file:
        reader = csv.DictReader(original_file, delimiter='\t')
        # write a header with the column names in the input file
        new_file.write('\t'.join(reader.fieldnames) + '\n')
        for line in reader:
            # skip entries associated with other species
            if line['Species'] in species:
                continue
            # skip entries without sufficient evidence
            if line['Evidence Subtypes'] in ['UNKNOWN', 'N/A']:
                continue
            new_file.write('\t'.join([line[c] for c in reader.fieldnames]) + '\n')


class Entry:
    """ MiRBase entry used for processing the miRBase data file. """
    def __init__(self, identifier: str, accession: str, species: str, prefix: str, ids: Set[str]):
        # store the basic information
        self.identifier = identifier            # type: str
        self.accession = accession              # type: str
        self.species = species                  # type: str
        # initialize the identifier set and add the accession number and the identifier to it, remove empty identifiers
        self.ids = ids                          # type: Set[str]
        self.ids.add(identifier)
        self.ids.add(accession)
        self.ids.discard('')
        # the main identifier
        self.symbol = identifier                # type: str
        # check which of the identifiers is actually the proper main identifier
        for x in ids:
            x = x.replace(prefix, '')
            if x.isdigit() or x[:3].lower() not in ['mir', 'let']:
                continue
            if x.replace('-', '').upper() in ids:
                self.symbol = x.replace('-', '').upper()
            elif x.lower().replace('mir-', 'mir').replace('let-', 'let').upper() in ids:
                self.symbol = x.lower().replace('mir-', 'mir').replace('let-', 'let').upper()


class MiRBaseParser:
    def __init__(self, data_path: str):
        """ Parses and processes a MiRBase data file to extract all identifier mappings and sort them by species. """
        self.data_path = data_path                          # type: str
        # key: species name (e.g. Homo sapiens), value: dictionary (key: main identifier, value: mapping entry)
        self.entries = defaultdict(dict)                    # type: Dict[str, Dict[str, Entry]]
        # extract all entries from the data file, process them and add them to the entry list
        with open(self.data_path, 'r') as file:
            entry_lines = []                                # type: List[str]
            for line in file:
                # end/start of a new entry => process the current entry, start an empty one
                if line.startswith('//'):
                    self.mirna_entry(entry_lines)
                    entry_lines = []
                # the line belongs to the current entry => add it to the entry
                else:
                    entry_lines.append(line.strip())

    def mirna_entry(self, lines: List[str]):
        """ Processes a list of lines belonging to a single entry in the MiRBase data file and stores it. """
        # nothing to do if the entry is empty
        if not lines:
            return
        # main information contained in the entry
        identifier = ''                                     # type: str
        accession = ''                                      # type: str
        species = ''                                        # type: str
        prefix = ''                                         # type: str
        ids = set()                                         # type: Set[str]
        # go over the entry line by line to extract the information
        for line in lines:
            tag, line = (line[:2], line[2:].strip())
            # identifier tag
            if tag == 'ID':
                identifier = line.split()[0]
                prefix = line.split(';')[2].strip().lower() + '-'
                ids.add(identifier.replace(prefix, ''))
                ids.add(identifier.replace(prefix, '').replace('mir-', 'miR-'))
                ids.add(identifier.replace('mir-', 'miR-'))
            # access tag
            elif tag == 'AC':
                accession = line.strip(';')
            # species tag
            elif tag == 'DE':
                species = line.split(identifier.replace(prefix, '').replace('mir', 'miR'))[0].strip()
            # further identifier(s) tag
            elif tag == 'DR':
                source, source_id, symbol = line.split('; ')
                symbol = symbol[:-1] if symbol.endswith('.') else symbol
                if source == 'HGNC':
                    ids.add('HGNC:' + source_id.strip())
                elif source == 'ENTREZGENE':
                    ids.add(source_id.strip())
                ids.add(symbol)
                ids.add(symbol.replace(prefix + 'miR-', prefix + 'mir-'))
        # the entry did not contain any valid information
        if not identifier:
            return
        # further process the information and consolidate them in an entry
        entry = Entry(identifier, accession, species, prefix, ids)
        # do not save the entry if there is some issue with the main identifier (symbol)
        if not entry.symbol:
            return
        # save the entry under its species and main identifier
        self.entries[species][entry.symbol] = entry

    def write_mapping(self, file_path: str, species: str):
        """
        Generates a .tsv file containing the MiRBase identifier mapping for the given species.
        :param file_path: output path for the species-specific MiRBase mapping
        :param species: scientific name (e.g. Homo sapiens) of the species of interest
        """
        # nothing to do if there is no mapping available for the species
        if species not in self.entries:
            return
        # write all mapping entries to the specified output file
        with open(file_path, 'w') as file:
            file.write('miRNA\tidentifiers\n')
            for key, entry in sorted(self.entries[species].items()):
                file.write('{0}\t{1}\n'.format(key, ', '.join(sorted(entry.ids))))


def gene_type_annotations(mirna_path: str, tf_path: str, mirbase_path: str, final_path: str):
    """
    This function compiles gene type annotations from
    :param mirna_path:
    :param tf_path:
    :param mirbase_path:
    :param final_path: path to the file that is going to contain the compiled information
    """
    print('Gene Type Annotation Compilation')
    with open(mirbase_path, 'r') as file:
        mirbase = {x for line in csv.DictReader(file, delimiter='\t')
                   for c in line.keys() if line[c].strip() for x in line[c].split(', ')}

    print('\tmirBase:', len(mirbase))

    with open(mirna_path, 'r') as file:
        mirna = {line[c] for line in csv.DictReader(file, delimiter='\t')
                 for c in line.keys() if c != 'Gene type' and line[c].strip()}

    print('\tEnsembl:', len(mirna))
    print('\tEnsembl mirBase overlap:', len(mirna.intersection(mirbase)))

    with open(tf_path, 'r') as file:
        tfs = {line[c] for line in csv.DictReader(file, delimiter=',')
               for c in ['Ensembl ID', 'HGNC symbol', 'EntrezGene ID'] if line[c].strip() and line['Is TF?'] == 'Yes'}

    # give some information about the number of gene types
    all_mirna = mirna.union(mirbase)
    print('\tTFs:', len(tfs))
    print('\tTF-miRNA overlap:', len(all_mirna.intersection(tfs)))
    print('\tTF-Ensembl overlap:', len(mirna.intersection(tfs)))
    print('\tTF-mirBase overlap:', len(mirbase.intersection(tfs)))

    # save the gene type annotations
    with open(final_path, 'w') as file:
        file.write('gene\tgene type\n')
        for m in sorted(all_mirna):
            file.write(f'{m}\tmiRNA\n')
        for t in sorted(tfs):
            if t not in all_mirna:
                file.write(f'{t}\tTF\n')

    print()


def compile_field_names(dir_path: str):
    """
    Compiles and prints the column names contained in all files in the raw mapping directory. Useful for setting up the
    mapping configuration in database.yaml
    :param dir_path: path to the directory with the raw mapping files
    """
    if not os.path.isdir(dir_path):
        print(f'FIELD NAMES: "{dir_path}" is not an existing directory')
        return

    names = set()
    for file_name in os.listdir(dir_path):
        with open(os.path.join(dir_path, file_name), 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter='\t')
            names.update(reader.fieldnames)

    names = sorted(names)
    print(f'FIELD NAMES: {names}')
    for x in sorted(names):
        print('* ' + x)
