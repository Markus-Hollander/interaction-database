from collections import defaultdict, Counter
import csv
from datetime import date
from itertools import combinations, takewhile
import sqlite3
from typing import Callable

from config import JobConfig
from database_components import Gene, Regulation, Symbol
from functions import *


class DataBaseBuilder:
    def __init__(self, cfg: JobConfig):
        # CONFIGURATION
        self.cfg = cfg                                                      # type: JobConfig

        # GENE SYMBOLS
        # names of symbols that were loaded from the main mapping
        self.symbols_main = set()                                           # type: set[str]
        # intermediary set of genes to collect gene specific information
        self.symbols_raw = defaultdict(Gene)                                # type: dict[str, Gene]
        # final set of genes with gene specific information, derived from the set of symbols after all mapping,
        # regulatory and annotation information has been processed and resolved
        self.genes = {}                                                     # type: dict[str, Gene]
        self.genes_discarded = dict()                                       # type: dict[str, Gene]

        # MAPPING
        # identifier: {symbol_name: Symbol object}
        self.mapping_raw = defaultdict(lambda: defaultdict(Symbol))         # type: dict[str, dict[str, Symbol]]
        # identifier: {sources}
        self.mapping_single = defaultdict(set)                              # type: dict[str, set[str]]
        #
        self.mapping_discarded = dict()
        # transformed identifier: {identifiers that were transformed into the key}
        self.transformed = defaultdict(set)                                 # type: dict[str, set[str]]
        self.ambiguous = {}
        self.mapping = {}

        # REGULATION
        # {(regulator, target): Regulation object}
        self.regulation_raw = defaultdict(Regulation)                       # type: dict[tuple[str, str], Regulation]
        self.regulation = defaultdict(Regulation)                           # type: dict[tuple[str, str], Regulation]
        self.regulation_single = defaultdict(set)                           # type: dict[str, set[str]]
        self.regulation_mapping = defaultdict(set)                          # type: dict[tuple[str, str], set[float]]
        self.regulation_type = defaultdict(set)                             # type: dict[str, set[str]]
        self.regulation_score = defaultdict(set)                            # type: dict[tuple[str, str], set[float]]
        self.ppi_only = set()                                               # type: set[str]

        # ANNOTATIONS
        # {annotation type: {annotation: {gene name 1, gene name 2,...}}}
        self.annotations_raw = defaultdict(lambda: defaultdict(set))        # type: dict[str, dict[str, set[str]]]
        self.annotations = defaultdict(lambda: defaultdict(set))            # type: dict[str, dict[str, set[str]]]

        # LOG
        self.log_dict = defaultdict(lambda: defaultdict(dict))              # type: dict[str, dict[str, dict]]

        # RUN
        if self.exec('mapping', self.load_mapping, has_sub=True):
            return
        if self.exec('regulation', self.load_regulation, flag=self.cfg.run.skip_regulation, has_sub=True):
            return
        if self.exec('resolve mapping', self.resolve_mapping):
            return
        if self.exec('load annotations', self.load_annotations, has_sub=True):
            return
        if self.exec('resolve regulation', self.resolve_regulation, flag=self.cfg.run.skip_regulation_resolve):
            return
        if self.exec('final cleanup', self.final_cleanup, flag=self.cfg.run.skip_final_cleanup):
            return
        if self.exec('build database', self.build_database, flag=self.cfg.run.skip_database, has_sub=True):
            return
        self.exec('report', self.report, flag=self.cfg.run.skip_log)

    # HELPER FUNCTIONS - RUN TIME INFORMATION --------------------------------------------------------------------------
    def error(self, indent: str, prefix: str, fmt: str) -> bool:
        """
        :return: True if a major error occurred, False otherwise
        """
        if self.cfg.run.skip:
            print(fmt.format(f'{indent}{prefix}ERROR:'), self.cfg.run.skip)
            return True
        return False

    def exec(self, label: str, func: Callable, flag='', has_sub=False, is_sub=False, timed=True, *args) -> bool:
        """
        :param label: message to print when calling the function
        :param func: function to execute
        :param flag: empty string if the function should be called, otherwise error message
        :param has_sub: True if the function call has sub-steps
        :param is_sub: True if the function is a sub-step of another function
        :param timed: True if the function should be timed, False otherwise
        :param args: function arguments
        :return: True if a major error occurred, False otherwise
        """
        # formatting
        indent = '        ' if is_sub else '    '           # type: str
        prefix = '' if is_sub else '--> '                   # type: str
        fmt = '{0:<42}'                                     # type: str
        # check for a major error
        if self.error(indent, label, fmt):
            return True
        # some run or error flag is set => do not execute the function
        if flag:
            print(fmt.format(f'{indent}{prefix}{label}:'), flag)
            return False
        # if the function has sub-steps, already print that the function is running
        if has_sub:
            print(f'{indent}{label}:')
        # execute the function
        start = time()
        return_value = func(*args)
        # print the finish message
        if timed:
            print(fmt.format(f'{indent}{prefix}{label}:'), f'{time_stamp(start)}')
        else:
            print(f'{indent}{prefix}{label}')
        # some issue occurred during a sub-step
        if return_value:
            return True
        # check once more for a major error
        return self.error(indent, label, fmt)

    # HELPER FUNCTIONS - LOGGING
    def log(self, category: str, sub_category: str, entry: str, key=None, update=None, add=None):
        """
        Adds an entry to the log dictionary, which is organised as follows:
        {category: {sub-category: {entry: value}}}
        The values can be lists, sets, dictionaries or numbers.
        :param category: dictionary key for the category
        :param sub_category: dictionary key for the sub-category
        :param entry: dictionary key for an entry in the sub-category
        :param key: key for the given value, only given if the entry is a dictionary
        :param update: given if the value is a string, set or list
        :param add: given if the value is a number
        """
        # make sure that exactly one value is given as input
        if update is None and add is None:
            raise ValueError('missing value')
        if update is not None and add is not None:
            raise ValueError('both update and add are given')
        # retrieve the entries of the specified sub-category
        entries = self.log_dict[category][sub_category]
        # if the key is given, the entry is a dictionary of sets
        if key:
            # if the entry does not exist yet, initialise it as a dictionary of sets
            if entry not in entries:
                entries[entry] = defaultdict(set)
            # from now on work on the value dictionary
            entries = entries[entry]
            entry = key
        # the value is a string, list or set that should get merged into the entry
        if update is not None:
            if entry not in entries:
                entries[entry] = set()
            if isinstance(update, str):
                entries[entry].add(update)
            else:
                entries[entry].update(update)
        # the value is a number that should get added to the entry
        elif add is not None:
            if entry in entries:
                entries[entry] += add
            else:
                entries[entry] = add

    def get(self, *keys, default=None):
        """
        Retrieve a value from the log dictionary.
        For example, get('a', 'b', 'c') will try to return log_dict['a']['b']['c'].
        :param keys: one or more keys
        :param default: value to return if the keys are not all in the log dictionary
        :return: specified log dictionary entry
        """
        # start at the top level of the log dictionary
        d = self.log_dict
        # iterate over the given keys to traverse the log dictionary until the last key is found, or until a key no
        # longer exists
        for i in range(len(keys)):
            # the key does not exist => return the default value
            if keys[i] not in d:
                return default
            d = d[keys[i]]
        # found the value
        return d

    # HELPER FUNCTIONS - IDENTIFIER PROCESSING -------------------------------------------------------------------------
    def find_existing_identifier(self, identifier: str) -> str:
        """
        Tries to find the identifier in various identifier collections, starting with the most refined one. For each
        collection, also try various spellings (upper case, lower case, capitalized) before looking in the next
        collection to avoid duplicate identifiers merely because of different cases.
        """
        for collection in [self.genes, self.mapping, self.symbols_main, self.symbols_raw, self.mapping_raw]:
            # continue with the next collection if the current one is empty
            if not collection:
                continue
            # check if the input identifier is known
            if identifier in collection:
                return identifier
            # if not, try different spellings
            if identifier.upper() in collection:
                return identifier.upper()
            if identifier.lower() in collection:
                return identifier.lower()
            if identifier.capitalize() in collection:
                return identifier.capitalize()
        # simply return the input identifier if it has not been added to any collection yet
        return identifier

    def process_identifier(self, identifier: str) -> tuple[str, set[str]]:
        """
        Checks the identifier for validity and performs some adjustments.
        :return: processed identifier, set of intermediate forms of the identifier
        """
        # remove excess whitespace
        identifier = identifier.strip()
        transformed_ids = {identifier}
        if identifier.lower() == 'gga-mir-9' or 'gg-mir-9' in identifier:
            print('AAAH: 0', identifier)
        # case: some notable error with the identifier
        if error := is_invalid(identifier):
            self.log('Identifier processing', 'invalid', error, update=identifier)
            return '', set()
        # case: GMT identifier
        if 'v$' in identifier.lower():
            transformed_ids.add(identifier)
            self.log('Identifier processing', 'invalid', 'GMT identifier', update=identifier)
            identifier = identifier.split('_')[1][2:]
        if '_' in identifier:
            transformed_ids.add(identifier)
            self.log('Identifier processing', 'adjustment', 'underscore(s)', update=identifier)
            identifier = identifier.replace('_', '-')
        if ';' in identifier:
            transformed_ids.add(identifier)
            self.log('Identifier processing', 'adjustment', 'semicolon(s)', update=identifier)
            identifier = identifier.replace(';', ',')
        # case: odd symbols in identifiers
        adjusted = ''.join([c for c in identifier if c.isalnum() or c in ['-', '+', ',', '(', ')']])
        if adjusted != identifier:
            transformed_ids.add(identifier)
            self.log('Identifier processing', 'adjustment', 'other odd character(s)', update=identifier)
            identifier = adjusted
        if identifier.startswith('(') and identifier.endswith(')'):
            transformed_ids.add(identifier)
            self.log('Identifier processing', 'adjustment', 'surrounding parentheses', update=identifier)
            identifier = identifier[1:-1]
        # case: various forms of miRNA identifiers
        if is_miRBase_accession(identifier):
            if identifier.lower() == 'gga-mir-9':
                print('AAAH: 1')
            self.log('Identifier processing', 'miRBase accession', 'raw', update=identifier)
            if identifier != identifier.upper():
                transformed_ids.add(identifier)
                self.log('Identifier processing', 'miRBase accession', 'new', update=identifier)
                identifier = identifier.upper()
        elif is_mirna_symbol(identifier):
            if identifier.lower() == 'gga-mir-9':
                print('AAAH: 2')
            self.log('Identifier processing', 'miRNA symbol', 'raw', update=identifier)
            if identifier != identifier.upper():
                transformed_ids.add(identifier)
                self.log('Identifier processing', 'miRNA symbol', 'new', update=identifier)
                identifier = identifier.upper()
        elif is_miRBase_ID(identifier, self.cfg.miRNA_prefix):
            self.log('Identifier processing', 'miRBase ID', 'raw', update=identifier)
            new_identifier, transformed = process_mirBase_ID(identifier, self.cfg.miRNA_prefix)
            if identifier.lower() == 'gga-mir-9':
                print('AAAH: 3', 'new id:', new_identifier, transformed)
            if not new_identifier:
                self.log('Identifier processing', 'invalid', 'wrong miRNA species', update=identifier)
                return '', set()
            transformed_ids.update(transformed)
            if identifier != new_identifier:
                self.log('Identifier processing', 'miRBase ID', 'new', update=identifier)
                identifier = new_identifier
        # check if the identifier is already known in a different case (upper, lower, capitalized)
        if identifier != (new_identifier := self.find_existing_identifier(identifier)):
            self.log('Identifier processing', 'found', 'raw', update=identifier)
            identifier = new_identifier
            self.log('Identifier processing', 'found', 'new', update=identifier)
        # remove the latest version of the identifier from the transformed IDs, in case it ended up in there
        transformed_ids.discard(identifier)

        return identifier, transformed_ids

    def add_raw_identifier(self, i: str, symbol: str, sources: str | Iterable[str], insensitive=False,
                           process_identifier=False, add_transformed=False, manual_fix=False, main=False,
                           main_fix=False, miRBase=False, entrez=False, regulation=False):
        """
        Adds an identifier - symbol pair to the raw mapping and performs some additional processing if specified.
        :param i: identifier
        :param symbol: gene symbol
        :param sources: single source or list of sources supporting the identifier - symbol mapping
        :param insensitive: True if the identifier spelling is case-insensitive, False otherwise
        :param process_identifier: True if the identifier should be processed prior to adding the mapping
        :param add_transformed: True if the intermediate versions of the identifier should be added as well
        :param manual_fix: True if the mapping is the result of a manual fix, False otherwise
        :param main: True if the mapping was part of the main mapping, False otherwise
        :param main_fix: True if the mapping was part of fixing the main mapping (untangle_main_mapping())
        :param miRBase: True if the mapping was taken from the miRBase mapping, False otherwise
        :param entrez: True if the mapping was taken from the Entrez mapping, False otherwise
        :param regulation: True if the mapping was added from regulatory information
        """
        # remove whitespace
        i = i.strip()
        transformed = set()
        # if specified, process the identifier before adding the mapping
        if process_identifier:
            i, transformed = self.process_identifier(i)
        # case: empty or invalid identifier
        if not i:
            return
        # add various capitalization of the identifier to the mapping with the symbol
        identifiers = {i, i.upper(), i.lower(), i.capitalize()} if insensitive else {i}
        # add the intermediate identifier versions to the mapping with the symbol as well
        if add_transformed:
            identifiers.update(transformed)
        # prepare the set of sources for the mapping
        sources = {sources.strip()} if isinstance(sources, str) else sources
        sources.discard('')
        # case: the identifier was uncategorized => add it to the list of single identifiers for later processing
        if symbol == self.cfg.no_symbol:
            for i in identifiers:
                self.mapping_single[i].update(sources)
            return
        # map the identifier(s) to the symbol
        for i in identifiers:
            self.mapping_raw[i][symbol].name = symbol
            self.mapping_raw[i][symbol].sources.update(sources)
            self.symbols_raw[symbol].name = symbol
            self.symbols_raw[symbol].aliases[i].update(sources)
            # flag the symbol (and with the added mapping) with the type of source backing it
            if main:
                self.mapping_raw[i][symbol].in_main = True
                self.symbols_main.add(symbol)
            if regulation:
                self.mapping_raw[i][symbol].in_regulation = True
            if manual_fix:
                self.mapping_raw[i][symbol].in_manual_fix = True
            if main_fix:
                self.mapping_raw[i][symbol].in_main_cleanup = True
            if entrez:
                self.mapping_raw[i][symbol].in_Entrez = True
            if miRBase:
                self.mapping_raw[i][symbol].in_miRBase = True

    def remove_raw_mapping(self, identifier: str, symbol: str):
        """ Remove an identifier - symbol pair from the raw mapping. """
        # make sure first, that the mapping exists
        if symbol not in self.symbols_raw.keys() or identifier not in self.mapping_raw.keys():
            self.log('Identifier processing', 'remove raw mapping', 'mapping does not exist',
                     update=(identifier, symbol))
            return
        # remove the identifier from the symbol aliases and the symbol from the identifier
        del self.symbols_raw[symbol].aliases[identifier]
        del self.mapping_raw[identifier][symbol]
        self.log('Identifier processing', 'remove raw mapping', 'removed mapping', update=(identifier, symbol))
        # if the symbol no longer has any aliases, remove it entirely
        if not self.symbols_raw[symbol].aliases:
            del self.symbols_raw[symbol]
            self.log('Identifier processing', 'remove raw mapping', 'removed symbol',
                     update=(identifier, symbol))
            # also remove it from the list of main symbols if it is in there
            if symbol in self.symbols_main:
                self.symbols_main.discard(symbol)
                self.log('Identifier processing', 'remove raw mapping', 'removed symbol, main',
                         update=(identifier, symbol))

    def process_aliases(self, input_aliases: Iterable[str], log_key: str, second_key: str, symbol='',
                        entrez='', allow_numerical=False) -> set[str]:
        """
        Processes and validates a set of input aliases.
        :param input_aliases: set of aliases that need to be processed and validated
        :param log_key: name of the log category
        :param second_key: name of the log sub-category
        :param symbol: gene symbol given alongside the alias column, empty string if the source does not contain one
        :param entrez: Entrez ID given alongside the alias column, empty string if the source does not contain one
        :param allow_numerical: True if aliases are allowed to be numerical (i.e. the alias column contains Entrez IDs
                                and not identifiers of other ID systems), False otherwise
        :return: set of processed and validated aliases
        """
        # make sure that the alias column is not empty
        if not input_aliases:
            self.log(log_key, second_key, 'aliases, empty column', add=1)
            return set()
        # otherwise process the alias(es)
        processed_aliases = set()
        for alias in sorted(input_aliases):
            # pre-process the identifier
            alias, transformed = self.process_identifier(alias)
            # invalid alias
            if not alias:
                self.log(log_key, second_key, 'aliases, invalid', add=1)
                continue
            # the alias is numerical and therefore potentially an Entrez ID
            if alias.isdigit():
                # log some information about the alias
                if entrez:
                    case = 'matching Entrez ID' if alias == entrez else 'not matching Entrez ID'
                    self.log(log_key, second_key, f'aliases, numerical ({case})', add=1)
                else:
                    self.log(log_key, second_key, f'aliases, numerical', add=1)
                # discard the alias if there should be no numerical aliases in the column (e.g. because the aliases
                # in that column could belong to a different numerical identifier system), or if the alias is the same
                # as the given Entrez ID (to avoid self-mapping)
                if not allow_numerical or alias == entrez:
                    continue
            # the alias in the symbol column is also contained in the alias column => skip the alias
            if symbol and alias == symbol:
                self.log(log_key, second_key, 'aliases, contains symbol', add=1)
                continue
            # the alias seems to be valid => add it to the list
            processed_aliases.add(alias)
            # only do this if there are actually transformed aliases, to avoid creating an entry in the default
            # dictionary if there are none
            if transformed:
                self.transformed[alias].update(transformed)
        return processed_aliases

    def add_additional_mapping(self, main_identifier: str, identifiers: set[str], input_symbols: set[str], source: str,
                               log_key: str, prefix: str, main=False, miRBase=False, entrez=False,
                               regulation=False) -> str | set[str]:
        """
        This function adds an identifier -> symbol mapping from additional sources to the main mapping from Ensembl.
        :param main_identifier:
        :param identifiers: set of identifiers associated with the main identifier (= aliases)
        :param input_symbols: set of gene symbols associated with the main identifier
        :param source: name of the mapping source
        :param log_key: name under which to log the processing information
        :param prefix: prefix for the logging messages
        :param main: True if the mapping source is the main mapping file, False otherwise
        :param miRBase: True if the mapping source is the miRBase mapping file, False otherwise
        :param entrez: True if the mapping source is the Entrez mapping file, False otherwise
        :param regulation: True if the mapping source is one of the regulation data files, False otherwise
        :return: set of valid and already known symbols contained in the input identifier set;
                 empty set => no known and valid symbols for the identifier
                 set of one symbol => unambiguous mapping
                 set of more than one symbol => ambiguous mapping
        """
        # make sure the main identifier is not in the set of associated identifiers
        identifiers.discard(main_identifier)
        # if there are no aliases for the main identifier, add the main identifier to the list of unmapped identifiers
        if not identifiers:
            self.log(log_key, f'{prefix} mapping', 'single', add=1)
            self.mapping_single[main_identifier].add(source)
            return set()

        # log if some symbols were given as input
        if not input_symbols:
            self.log(log_key, f'{prefix} ID processing', 'IDs from symbol column, none', add=1)
        elif len(input_symbols) == 1:
            self.log(log_key, f'{prefix} ID processing', 'IDs from symbol column, single', add=1)
        else:
            self.log(log_key, f'{prefix} ID processing', 'IDs from symbol column, multiple', add=1)

        # for further processing add the main identifier to the set of identifiers
        identifiers.add(main_identifier)

        # some logging regarding how many of the identifiers are already contained in the main mapping
        known_IDs = {i for i in identifiers if i in self.mapping_raw.keys()}
        if not known_IDs:
            self.log(log_key, f'{prefix} ID processing', 'known IDs, none', add=1)
        elif len(known_IDs) == len(identifiers):
            self.log(log_key, f'{prefix} ID processing', 'known IDs, all', add=1)
        else:
            self.log(log_key, f'{prefix} ID processing', 'known IDs, some', add=1)
        if known_IDs and main_identifier in known_IDs:
            self.log(log_key, f'{prefix} ID processing', 'known IDs, including main identifier', add=1)

        # get the subset of identifiers that are known gene symbols from the main mapping file
        known_symbols = {i for i in identifiers if i in self.symbols_main}
        if known_symbols:
            self.log(log_key, f'{prefix} ID processing', 'known symbols, directly from main mapping', add=1)
        # if none of the identifiers are known gene symbols from the main mapping file, check if any of the identifiers
        # are associated with symbols from the main mapping file
        else:
            known_identifiers = {i for i in identifiers if i in self.mapping_raw.keys()}
            known_symbols = {s for i in known_identifiers for s in self.mapping_raw[i].keys() if s in self.symbols_main}
            if known_symbols:
                self.log(log_key, f'{prefix} ID processing', 'known symbols, indirectly from main mapping',
                         add=1)
            else:
                # if that is not the case, continue with all symbols associated with the known identifiers
                known_symbols = {s for i in known_identifiers for s in self.mapping_raw[i].keys()}
                if known_symbols:
                    self.log(log_key, f'{prefix} ID processing', 'known symbols, other source', add=1)
                else:
                    self.log(log_key, f'{prefix} ID processing', 'known symbols, none', add=1)

        # log information about the known symbols
        if known_symbols and main_identifier in known_symbols:
            self.log(log_key, f'{prefix} ID processing', 'known symbols, including main identifier', add=1)
        if known_symbols and input_symbols and known_symbols.intersection(input_symbols):
            self.log(log_key, f'{prefix} ID processing', 'known symbols, including symbol column', add=1)
        if known_IDs and not known_symbols:
            self.log(log_key, f'{prefix} ID processing', 'known IDs, but no known symbols', add=1)

        # see if there is overlap between the input symbols and the known symbols
        symbols = known_symbols.intersection(input_symbols)
        if symbols:
            self.log(log_key, f'{prefix} ID processing', 'symbols, used overlap input and known', add=1)
        # if not, continue with the known symbols
        elif known_symbols:
            self.log(log_key, f'{prefix} ID processing', 'symbols, used known', add=1)
            symbols = known_symbols
        # if there is no overlap and no known symbols, then take the input symbols if given
        elif input_symbols:
            self.log(log_key, f'{prefix} ID processing', 'symbols, used input', add=1)
            symbols = input_symbols
        # otherwise use the main identifier as the symbol
        else:
            self.log(log_key, f'{prefix} ID processing', 'symbols, used main identifier', add=1)
            symbols = {main_identifier}

        # case: unambiguous symbol => map all identifiers to that symbol
        if len(symbols) == 1:
            symbol = list(symbols)[0]
            for identifier in identifiers:
                # avoid self-mapping for now
                if identifier == symbol:
                    continue
                self.add_raw_identifier(identifier, symbol, source, main=main, entrez=entrez, miRBase=miRBase,
                                        regulation=regulation)
            self.log(log_key, f'{prefix} mapping', 'unambiguous', add=1)
            return {symbol}

        # case: ambiguous symbols => return the symbol candidates
        self.log(log_key, f'{prefix} mapping', 'ambiguous', add=1)

        return symbols

    # MAPPING ----------------------------------------------------------------------------------------------------------
    def load_mapping(self):
        """
        Process all mapping files that are enabled in the configuration.
        :return: True if an error occurred in a mapping step, False otherwise
        """
        # generate the mapping first by merging (Ensembl Biomart) files in the given directory
        if self.exec('merge mapping', self.generate_mapping, flag=self.cfg.run.skip_merge, is_sub=True):
            return True
        # start with manual fixes, e.g. to avoid automatic Excel conversion issues
        if self.exec('load manual fixes', self.load_fixes, flag=self.cfg.run.skip_fixes, is_sub=True):
            return True
        # load the main mapping (likely from Ensembl)
        if self.exec('load main mapping', self.load_main_mapping, is_sub=True):
            return True
        # conduct some cleanup of the main mapping
        if self.exec('untangle main mapping', self.untangle_main_mapping, flag=self.cfg.run.skip_untangle, is_sub=True):
            return True
        # load additional mapping from miRBase
        if self.exec('load miRBase mapping', self.load_mirbase_mapping, flag=self.cfg.run.skip_miRBase, is_sub=True):
            return True
        # load additional mapping from Entrez
        if self.exec('load Entrez mapping', self.load_additional_entrez, flag=self.cfg.run.skip_Entrez, is_sub=True):
            return True
        return False

    def generate_mapping(self, symbol='Gene name'):
        """
        Generate a single mapping file from a directory full of smaller mapping files. It is assumed that each column
        in the mapping files contains only one identifier per cell, and that each file contains one main column
        to which to map the identifiers in the other columns.
        :param symbol: name of the column that contains the gene symbol
        """
        # {gene symbol: {source: {id1, id2,...}}}
        entries = defaultdict(lambda: defaultdict(set))         # type: dict[str, dict[str, set[str]]]
        column_names = set()                                    # type: set[str]

        # process all files in the merge directory
        for file_name in sorted(os.listdir(self.cfg.paths.mapping_merge_dir)):
            with open(os.path.join(self.cfg.paths.mapping_merge_dir, file_name), 'r', encoding='utf-8') as file:
                reader = csv.DictReader(file, delimiter='\t')

                for line in reader:
                    # label the identifiers of rows without symbol column as uncategorized
                    if symbol not in line or not line[symbol]:
                        line[symbol] = self.cfg.no_symbol
                    # add the identifiers of each column (including the symbol column)
                    for col in reader.fieldnames:
                        # skip columns without identifier
                        if not line[col]:
                            continue
                        # store the symbol - column mapping and the column name
                        entries[line[symbol]][col].add(line[col])
                        column_names.add(col)
        # prepare the column order for the mapping file
        columns = [symbol]
        columns = columns + sorted(x for x in column_names if x not in columns)
        # output the mapping file
        with open(self.cfg.paths.main.path, 'w', encoding='utf-8') as file:
            # write the header
            file.write('\t'.join(columns) + '\n')
            # write the entries sorted by gene symbol
            for _, entry in sorted(entries.items()):
                cols = [', '.join(sorted(entry[col])) for col in columns]
                file.write('\t'.join(cols) + '\n')

    def load_fixes(self, symbol='gene.symbol', to_fix_col='to.fix', case='case.insensitive'):
        """
        Due to automatic data conversion in some data and table processing software, e.g. Excel, some data sets
        contain erroneous identifiers. This allows to load manual fixes to these and other errors. This should be done
        before loading the main mapping, in case the main mapping contains some of these errors.
        :param symbol: main column name, contains the gene symbol
        :param to_fix_col: column name, contains the identifiers that need to be converted to the main symbol
        :param case: column name, contains 'yes' if the 'to_fix' identifiers are case sensitive, 'no' otherwise
        """
        with open(self.cfg.paths.fixes.path, 'r') as file:
            for line in csv.DictReader(file, delimiter='\t'):
                # skip rows without gene symbol information
                if symbol not in line or not line[symbol]:
                    continue
                # map the symbol to itself
                self.add_raw_identifier(line[symbol], line[symbol], 'manual fix', manual_fix=True)
                # iterate over all identifiers that should get converted to the symbol
                for to_fix in line[to_fix_col].split(', '):
                    self.add_raw_identifier(to_fix, line[symbol], 'manual fix', insensitive=line[case] == 'yes',
                                            manual_fix=True)

        self.log('Initial mapping', 'manual fixes', 'symbols', add=len(self.symbols_main))
        self.log('Initial mapping', 'manual fixes', 'identifiers', add=len(self.mapping_raw))

    def load_main_mapping(self, symbol_col='Gene name'):
        """
        Read the identifier - symbol mapping definition file (likely from Ensembl).
        :param symbol_col: column with the gene symbol
        """
        # make sure the symbol is mapped to itself
        self.cfg.run.main_cols.add(symbol_col)
        # preloading stats
        symbols_main = len(self.symbols_main)
        identifiers_raw = len(self.mapping_raw)
        # load the mapping
        with open(self.cfg.paths.main.path, 'r', encoding='utf-8') as file:
            for line in csv.DictReader(file, delimiter='\t'):
                # flag entries without declared gene symbol
                if not line[symbol_col]:
                    self.log('Initial mapping', 'main', 'symbols, empty', add=1)
                # get the symbol
                symbol = line[symbol_col].strip()
                # for the other defined columns, add the mapping to the gene symbol
                for col in self.cfg.run.main_cols:
                    # skip the column if it is not in the file or if the cell is empty
                    if col not in line or not line[col]:
                        continue
                    # skip the column if it is the symbol column, and it is the catch-all symbol
                    if col == symbol_col and line[col] == self.cfg.no_symbol:
                        continue
                    # get all identifiers in the cell (comma separated)
                    for identifier in line[col].split(', '):
                        processed, transformed = self.process_identifier(identifier)
                        # case: invalid identifier
                        if not processed:
                            self.log('Initial mapping', 'main', 'identifiers, invalid', update=identifier)
                            continue
                        # case: numerical identifier in a column that is not marked as containing Entrez identifiers
                        # => skip to avoid confusion with other numerical ID systems
                        if processed.isdigit() and col not in self.cfg.run.Entrez_cols:
                            self.log('Initial mapping', 'main', 'identifiers, non-Entrez number', update=processed)
                            continue
                        # store the transformed identifiers for later processing
                        if transformed:
                            self.log('Initial mapping', 'main', 'identifiers, transformed', add=1)
                            self.transformed[processed].update(transformed)
                        # add the identifier - symbol mapping, which already takes care of the case: symbol = NO_SYMBOL
                        self.add_raw_identifier(processed, symbol, col, main=True)
        # log the number of added symbols and identifiers
        self.log('Initial mapping', 'main', 'symbols', add=len(self.symbols_main) - symbols_main)
        self.log('Initial mapping', 'main', 'identifiers, all added', add=len(self.mapping_raw) - identifiers_raw)

    def untangle_main_mapping(self):
        """
        Perform some initial untangling of the main mapping to make later additions less bothersome.
        """
        # MIRBASE IDs as SYMBOLS
        # miRBase_ID (e.g. hsa-mir-103) should ideally not be used as symbols
        # => iterate over the miRBase IDs in the symbol list
        for miRBase_ID in sorted({symbol for symbol in self.symbols_main
                                  if is_miRBase_ID(symbol, self.cfg.miRNA_prefix)}):
            self.log('Initial mapping', 'miRBase IDs as symbols', 'all', add=1)
            # the miRBase ID is also in the mapping
            if miRBase_ID in self.mapping_raw.keys():
                # get the symbols associated with the miRBase ID
                other_miRBase_IDs = set(self.mapping_raw[miRBase_ID].keys())
                # convert the miRBase ID to miRNA symbol and check if it is among the associated symbols
                overlap = other_miRBase_IDs.intersection(miRBase_ID_to_symbol(miRBase_ID, self.cfg.miRNA_prefix).keys())
                # if there is a unique match, map the miRBase ID to the match, and to the same for all aliases of
                # the miRBase ID, then remove the previous mappings
                if len(overlap) == 1:
                    ids = list(self.symbols_raw[miRBase_ID].aliases.keys())
                    for i in ids:
                        self.add_raw_identifier(i, list(overlap)[0], self.mapping_raw[i][miRBase_ID].sources,
                                                main_fix=True, main=True)
                        self.remove_raw_mapping(i, miRBase_ID)
                    self.log('Initial mapping', 'miRBase IDs as symbols', 'fixed', add=1)
                else:
                    self.log('Initial mapping', 'miRBase IDs as symbols', 'ambiguous', update=miRBase_ID)
            else:
                self.log('Initial mapping', 'miRBase IDs as symbols', 'not mapped to other', update=miRBase_ID)

        # AMBIGUOUS MAPPINGS
        # compile all identifier that map to more than one symbol and try to resolve the ambiguity
        ambiguous_ids = {identifier for identifier, symbols in self.mapping_raw.items() if len(symbols) > 1}
        self.log('Initial mapping', 'main', 'identifiers, ambiguous', add=len(ambiguous_ids))

        # key: tuple of symbols, values: set of identifiers that are ambiguously mapped to that combination of symbols
        accession_to_resolve = defaultdict(set)                                     # type: dict[tuple[str], set[str]]
        # process all ambiguously mapped identifiers
        for identifier in sorted(ambiguous_ids):
            # get all symbols associated with the identifier
            symbols = self.mapping_raw[identifier]
            # get all sources backing the mapping of the identifier to the symbols and count the occurrence of each
            # source
            sources = Counter([source for ids in symbols.values() for source in ids.sources])
            # case: miRBase ID mapping with ambiguous symbols, check if any of them is the miRNA symbol
            if 'miRBase ID' in sources:
                # convert the miRBase ID to miRNA symbol and check if it is among the associated symbols
                mirna_symbols = miRBase_ID_to_symbol(identifier, self.cfg.miRNA_prefix)
                overlap = list(set(symbols.keys()).intersection(mirna_symbols.keys()))
                self.log('Initial mapping', 'untangle miRBase IDs', 'ambiguous, all', add=1)
                # if there is a unique match, it is the miRNA symbol => choose that as the correct mapping for the
                # miRBase ID, then remove the mappings to the other symbols
                if len(overlap) == 1:
                    symbol = overlap[0]
                    self.mapping_raw[identifier][symbol].in_main_cleanup = True
                    for k in set(symbols.keys()).difference({symbol}):
                        self.remove_raw_mapping(identifier, k)
                    self.log('Initial mapping', 'untangle miRBase IDs', f'unambiguous, {mirna_symbols[symbol]}',
                             update=identifier)
                elif not overlap:
                    self.log('Initial mapping', 'untangle miRBase IDs', 'overlap, none', update=identifier)
                else:
                    self.log('Initial mapping', 'untangle miRBase IDs', 'overlap, multiple', update=identifier)
            # case: miRBase accessions with ambiguous symbols, save the identifier grouped by the symbol
            # combination for later processing
            elif 'miRBase accession' in sources:
                self.log('Initial mapping', 'untangle miRBase accessions', 'ambiguous, all', add=1)
                # if all symbols are miRBase accessions, the situation could maybe be resolved => save it
                if sources['miRBase accession'] == len(symbols):
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'ambiguous, all sources include "miRBase accession"', update=identifier)
                    accession_to_resolve[tuple(sorted(symbols.keys()))].add(identifier)
                else:
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'ambiguous, not all sources include "miRBase accession"', update=identifier)

        # AMBIGUOUS MIRBASE ACCESSIONS
        # this resolves ambiguity when a set if miRBase accessions and miRNA symbols are all mapped to each other.
        # for example: symbols = [MIR1, MIR2] and accessions = [MI0000001, MI0000002]
        for symbols, accessions in sorted(accession_to_resolve.items()):
            # if the number of accessions and symbols is not the same, remove non-miRBase accessions and symbols from
            # the respective sets to see if that resolves the unbalance
            if len(accessions) != len(symbols):
                # remove miRBase IDs from the symbols
                n = len(symbols)
                symbols = sorted({s for s in symbols if not is_miRBase_ID(s, self.cfg.miRNA_prefix)})
                if n > len(symbols):
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'all sources include "miRBase accession", remove miRBase IDs from symbols', add=1)
                # remove non-miRBase accession from the accessions
                n = len(accessions)
                accessions = {i for i in accessions if is_miRBase_accession(i)}
                if n > len(accessions):
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'all sources include "miRBase accession", remove non-miRBase accessions from accessions ',
                             add=1)
            # case: the number of miRBase accessions and miRNA symbols are mapped to each other is the same.
            # => can be untangled by matching the symbols and accessions in alphabetical order
            if len(accessions) == len(symbols):
                for symbol, accession in zip(symbols, sorted(accessions)):
                    self.mapping_raw[accession][symbol].in_main_cleanup = True
                    # remove the mapping of the miRBase accession to the other miRNA symbols to resolve the ambiguity
                    for k in set(self.mapping_raw[accession].keys()).difference({symbol}):
                        self.remove_raw_mapping(accession, k)
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'all sources include "miRBase accession", all tangled', update=accession)
                continue
            # case: a single miRBase accession is mapped to multiple miRNA symbols
            if len(symbols) > 1 and len(accessions) == 1:
                # get the miRBase accession
                accession = list(accessions)[0]
                # get the symbols that are not yet mapped to a miRBase accession (other than the miRBase accession
                # under consideration)
                not_yet = {s for s in symbols
                           if s in self.symbols_raw.keys()
                           and self.symbols_raw[s].get_miRBase_accessions().difference({accession})}
                # if there is only a single such miRNA symbol, map the miRBase accession to that symbol
                if len(not_yet) == 1:
                    symbol = list(not_yet)[0]
                    # remove the mapping of the miRBase accession to the other miRNA symbols to resolve the ambiguity
                    for k in set(self.mapping_raw[accession].keys()).difference({symbol}):
                        self.remove_raw_mapping(accession, k)
                    self.log('Initial mapping', 'untangle miRBase accessions',
                             'all sources include "miRBase accession", single unassociated', update=accession)
                    continue
            # if the ambiguity could not be resolved, log the accessions (without duplicates)
            for accession in accessions:
                self.log('Initial mapping', 'untangle miRBase accessions',
                         'all sources include "miRBase accession", unresolved', update=accession)

    def load_mirbase_mapping(self, mirna_col='miRNA', alias_col='identifiers'):
        """
        In addition to the main mapping from Ensembl, this loads miRNA centric mapping information from miRBase and
        integrates it with the main mapping.
        :param mirna_col: name of the column containing the primary miRBase identifier
        :param alias_col: name of the column containing a list of aliases associated with the miRNA symbol
        """
        # note down the current numbers to later determine and log the changes
        symbols_raw = len(self.symbols_raw)
        identifiers_raw = len(self.mapping_raw)
        mirbase_ids_as_symbols = len({s for s in self.symbols_raw.keys() if is_miRBase_ID(s, self.cfg.miRNA_prefix)})
        old_ambiguous_ids = len({k for k, v in self.mapping_raw.items() if len(v) > 1})
        # set of miRBase identifiers and corresponding alias mapping contained in the file
        mirbase_ids = set()                                                             # type: set[str]
        mapping = defaultdict(set)                                                      # type: dict[str, set[str]]

        # read and process the miRBase file
        with open(self.cfg.paths.miRBase.path, 'r') as file:
            for line in csv.DictReader(file, delimiter='\t'):
                self.log('Initial mapping', 'miRBase load', 'rows, all', add=1)
                # skip rows with empty columns
                if not line[mirna_col]:
                    self.log('Initial mapping', 'miRBase load', 'rows, empty identifier column', add=1)
                    continue
                miRBase = line[mirna_col].strip()
                # skip numerical identifiers, since those cannot be confirmed to be or be distinguished from Entrez IDs
                if miRBase.isdigit():
                    self.log('Initial mapping', 'miRBase load', 'rows, numerical identifier', add=1)
                    continue
                # process and validate the identifier, and skip it if it's invalid
                miRBase, transformed = self.process_identifier(miRBase)
                if not miRBase:
                    self.log('Initial mapping', 'miRBase load', 'rows, invalid identifier', add=1)
                    continue
                self.log('Initial mapping', 'miRBase load', 'rows, valid', add=1)
                # add the identifier to the list
                mirbase_ids.add(miRBase)
                if transformed:
                    self.log('Initial mapping', 'miRBase load', 'rows, valid - with transformed', add=1)
                    self.transformed[miRBase].update(transformed)
                # if aliases are given for the identifier, add them to the index for further processing
                if aliases := self.process_aliases(line[alias_col].split(','), 'Initial mapping', 'miRBase load'):
                    mapping[miRBase].update(aliases)

        # process the loaded mapping
        for miRBase in sorted(mirbase_ids):
            # process and add the new mapping to the main mapping
            symbols = self.add_additional_mapping(miRBase, mapping[miRBase], set(), 'miRBase', 'Initial mapping',
                                                  'miRBase', miRBase=True)
            # cases 1 and 2: only miRBase ID (no mapping required) or unambiguous symbols => already processed above
            if not symbols or len(symbols) == 1:
                continue
            # case 3: ambiguous symbols => still needs to be processed in some way
            if is_miRBase_ID(miRBase, self.cfg.miRNA_prefix):
                # check if there is a single miRNA symbol amongst the mapping and map the identifier to that symbol
                miRBase_symbols = miRBase_ID_to_symbol(miRBase, self.cfg.miRNA_prefix)
                overlap = list(symbols.intersection(miRBase_symbols.keys()))
                if not overlap:
                    self.log('Initial mapping', 'miRBase mapping', 'ambiguous, miRNA symbol - no match', add=1)
                elif len(overlap) == 1:
                    self.log('Initial mapping', 'miRBase mapping', 'ambiguous, miRNA symbol - single match', add=1)
                    self.add_raw_identifier(miRBase, overlap[0], 'miRBase', main=True)
                else:
                    self.log('Initial mapping', 'miRBase mapping', 'ambiguous, miRNA symbol - multiple matches', add=1)
            else:
                self.log('Initial mapping', 'miRBase mapping', 'ambiguous, not miRNA symbol', add=1)

        # log the mapping changes
        self.log('Initial mapping', 'miRBase mapping', 'miRBase IDs, total', add=len(mirbase_ids))
        self.log('Initial mapping', 'miRBase summary', 'new symbols, total', add=len(self.symbols_raw) - symbols_raw)
        self.log('Initial mapping', 'miRBase summary', 'identifiers, new',
                 add=len(self.mapping_raw) - identifiers_raw)
        self.log('Initial mapping', 'miRBase summary', 'new symbols, miRBase IDs', add=len(
            {s for s in self.symbols_raw.keys() if is_miRBase_ID(s, self.cfg.miRNA_prefix)}) - mirbase_ids_as_symbols)
        self.log('Initial mapping', 'miRBase summary', 'identifiers, new ambiguity',
                 add=len({k for k, v in self.mapping_raw.items() if len(v) > 1}) - old_ambiguous_ids)

    def load_additional_entrez(self, symbol_col='Gene name', entrez_col='Entrez ID', aliases_col='Entrez Synonyms'):
        """
        In addition to the main mapping from Ensembl, this loads  mapping information from NCBI Entrez and integrates it
        with the main mapping.
        :param symbol_col: name of the column containing the gene symbol
        :param entrez_col: name of the column containing the Entrez identifier
        :param aliases_col: name of the column containing a list of aliases associated with the Entrez identifier
        """
        # note down the current numbers to later determine and log the changes
        old_symbols_raw = len(self.symbols_raw)
        old_identifiers_raw = len(self.mapping_raw)
        old_ambiguous_ids = len({k for k, v in self.mapping_raw.items() if len(v) > 1})
        # set of Entrez identifiers, corresponding gene symbol and alias mapping contained in the file
        entrez_ids = set()                                                                  # type: set[str]
        entrez_symbols = defaultdict(set)                                                   # type: dict[str, set[str]]
        entrez_mapping = defaultdict(set)                                                   # type: dict[str, set[str]]

        # read and process the Entrez file
        with open(self.cfg.paths.Entrez.path, 'r') as file:
            for line in csv.DictReader(file, delimiter='\t'):
                # case: Entrez ID invalid (missing, not numerical)
                if not (entrez := line[entrez_col].strip()):
                    self.log('Initial mapping', 'Entrez load', 'Entrez, invalid', add=1)
                    continue
                if not entrez.isdigit():
                    self.log('Initial mapping', 'Entrez load', 'Entrez, not numerical', add=1)
                    continue
                # add the Entrez ID
                entrez_ids.add(entrez)
                # process the symbol column
                if symbol := line[symbol_col].strip():
                    symbol, transformed = self.process_identifier(symbol)
                    if symbol and symbol.isdigit():
                        case = 'matching Entrez ID' if symbol == entrez else 'not matching Entrez ID'
                        self.log('Initial mapping', 'Entrez load', f'symbol, numerical ({case})', add=1)
                    elif symbol:
                        entrez_symbols[entrez].add(symbol)
                        entrez_mapping[entrez].add(symbol)
                        if transformed:
                            self.transformed[symbol].update(transformed)
                    else:
                        self.log('Initial mapping', 'Entrez load', 'symbol, invalid', add=1)
                else:
                    self.log('Initial mapping', 'Entrez load', 'symbol, empty column', add=1)
                # process the column with other identifiers or accessions
                if identifiers := self.process_aliases(line[aliases_col].strip().split(', '),
                                                       'Initial mapping', 'Entrez load',
                                                       symbol, entrez, False):
                    entrez_mapping[entrez].update(identifiers)

        # process the loaded mapping
        for entrez in sorted(entrez_ids):
            symbols = self.add_additional_mapping(entrez, entrez_mapping[entrez], entrez_symbols[entrez],
                                                  'Entrez', 'Initial mapping', 'Entrez', entrez=True)
            # cases 1 and 2: only Entrez ID (no mapping required) or unambiguous symbols => already processed
            if not symbols or len(symbols) == 1:
                continue
            # case 3: ambiguous symbols => still needs to be processed in some way
            # case 3a: Entrez ID is not yet known => no way to disambiguate the symbols
            if entrez not in self.mapping_raw.keys():
                self.log('Initial mapping', 'Entrez mapping', 'ambiguous, Entrez ID unknown', add=1)
                continue
            # case 3b:
            overlap = symbols.intersection(self.mapping_raw[entrez].keys())
            if not overlap:
                self.log('Initial mapping', 'Entrez mapping', 'ambiguous, no Entrez ID overlap', add=1)
            elif len(overlap) == 1:
                self.log('Initial mapping', 'Entrez mapping', 'ambiguous, single Entrez ID overlap', add=1)
            else:
                self.log('Initial mapping', 'Entrez mapping', 'ambiguous, multiple Entrez ID overlap', add=1)

        # log the mapping changes
        self.log('Initial mapping', 'Entrez mapping', 'Entrez IDs, total', add=len(entrez_ids))
        self.log('Initial mapping', 'Entrez summary', 'new identifiers',
                 add=len(self.mapping_raw) - old_identifiers_raw)
        self.log('Initial mapping', 'Entrez summary', 'new symbols',
                 add=len(self.symbols_raw) - old_symbols_raw)
        self.log('Initial mapping', 'Entrez summary', 'identifiers, new ambiguity',
                 add=len({k for k, v in self.mapping_raw.items() if len(v) > 1}) - old_ambiguous_ids)

    # REGULATION -------------------------------------------------------------------------------------------------------
    def load_regulation(self):
        # extract regulation information from all specified files
        def extract():
            for r in self.cfg.paths.regulation:
                self.extract_regulation(source=r.source, evidence=r.evidence, file_path=r.path,
                                        r_cols=r.regulator_cols, r_type=r.regulator_type, t_cols=r.target_cols,
                                        t_type=r.target_type, delimiter=r.delimiter, header=r.header,
                                        remaining_cols=r.remaining_cols, score_col=r.score_col,
                                        score_divisor=r.score_divisor, species_cols=r.species_cols, ppi_only=r.ppi_only)
        if self.exec('extract regulation', extract, is_sub=True):
            return True
        # process and untangle the mapping extracted from the loaded regulation data
        if self.exec('process regulation mapping', self.process_regulation_mapping,
                     flag=self.cfg.run.skip_regulation_untangle, is_sub=True):
            return True
        # note if no regulatory information could be loaded, which means the final regulation processing will be
        # skipped later
        if not self.cfg.run.skip_regulation_resolve and not self.regulation_raw:
            self.cfg.run.skip_regulation_resolve = 'no regulatory information'

        return False

    def get_identifier(self, source: str, first: str, second=''):
        """
        Takes a pair of mapped identifiers from regulation data and tries to determine which one is the primary one
        (i.e. most likely to be a gene symbol).
        :param source: name of the source
        :param first: first identifier mentioned in the source
        :param second: second identifier mentioned in the source
        :return: primary identifier
        """
        # PRE-PROCESSING
        first, first_transformed = self.process_identifier(first)
        second, second_transformed = self.process_identifier(second)

        # both identifiers are empty => cannot map anything
        if not first and not second:
            self.log('Regulation initial', 'identifier extraction', 'both empty', add=1)
            return ''

        # store the transformed IDs
        if first_transformed:
            self.transformed[first].update(first_transformed)
            self.log('Regulation initial', 'identifier extraction', 'first has transformed', add=1)
        if second_transformed:
            self.transformed[second].update(second_transformed)
            self.log('Regulation initial', 'identifier extraction', 'second has transformed', add=1)

        # CASES: only single identifier
        # both identifiers are the same => same as if only one identifier is given
        if first == second:
            self.log('Regulation initial', 'identifier extraction', 'same', add=1)
            second = ''

        # the secondary identifier is available but no the first  one => switch around
        if not first and second:
            self.log('Regulation initial', 'identifier extraction', 'first not given', add=1)
            first, second = (second, first)

        # only one of the identifiers is given => cannot establish a new mapping
        if not second:
            # if the identifier is unknown, save it for later integration
            if first not in self.mapping_raw.keys():
                if first not in self.mapping_single:
                    self.regulation_single[first].add(source)
                else:
                    self.log('Regulation initial', 'identifier extraction', 'single, already in single mapping', add=1)
            else:
                self.log('Regulation initial', 'identifier extraction', 'single, known', add=1)
            return first

        # both are valid => store the pair for later integration into the mapping
        self.regulation_mapping[tuple(sorted([first, second]))].add(source)

        return first

    def extract_regulation(self, source: str, evidence: str, file_path: str, r_cols, t_cols, delimiter='\t',
                           header=True, remaining_cols=False, r_type=None, t_type=None, score_col=None,
                           score_divisor=1.0, species_cols=None, ppi_only=False):
        """
        Extracts regulatory (and potentially additional mapping) information from a file.
        :param source: name of the regulation source
        :param evidence: type of the evidence (i.e. experimental or predicted)
        :param file_path: path to the file containing the regulatory information
        :param r_cols: list of column names (or indices) containing the identifier(s) of the regulator
        :param t_cols: list of column names (or indices) containing the identifier(s) of the target gene
        :param delimiter: delimiter between the columns
        :param header: True if the file contains a header row (with column names), False otherwise
        :param remaining_cols: True if the target gene(s) are given in the columns after the given target column index
        :param r_type: regulator type (i.e. TF, miRNA, gene)
        :param t_type: target type (i.e. TF, miRNA, gene)
        :param score_col: name (or index) of the column containing the interaction score
        :param score_divisor: number by which to divide the score
        :param species_cols: list of columns that contain species information for the regulation
        :param ppi_only: True if the regulation from this file is on the PPI level only
        """
        # keep track of PPI only sources
        if ppi_only:
            self.ppi_only.add(source)
        # consider the first column as the primary one
        r_col, r_col_secondary = (r_cols[0], r_cols[1:]) if len(r_cols) > 1 else (r_cols[0], [])
        t_col, t_col_secondary = (t_cols[0], t_cols[1:]) if len(t_cols) > 1 else (t_cols[0], [])
        # species to check
        species_to_check = [self.cfg.species.lower(), self.cfg.species_full.lower(), self.cfg.miRNA_prefix.lower()]
        # process the file
        with open(file_path, 'r') as file:
            # if the header is given, read each line as a dictionary to access the columns by name, otherwise read it
            # line by line and split it manually with the delimiter to access the column by index
            reader = csv.DictReader(file, delimiter=delimiter) if header else file
            # process each line
            for line in reader:
                # either a list of columns or a dictionary of columns
                line = line if header else line.strip().split(delimiter)

                # SPECIES CHECK
                if species_cols:
                    wrong_species = False
                    # check each column with species information if the species is consistent with the current species
                    for col in species_cols:
                        if line[col].lower() not in species_to_check:
                            wrong_species = True
                    # if there is any mismatch, ignore this entry
                    if wrong_species:
                        self.log('Regulation initial', 'loading', 'wrong species', add=1)
                        continue

                # REGULATOR EXTRACTION
                # if more than one regulator column is given, process the proposed mapping and determine the primary
                # identifier of each pair
                if r_col_secondary:
                    regulators = {self.get_identifier(source, line[r_col], line[secondary])
                                  for secondary in r_col_secondary}
                else:
                    regulators = {self.get_identifier(source, line[r_col])}
                # make sure at least one valid regulator identifier is included in the line
                if not regulators:
                    self.log('Regulation initial', 'loading', 'regulator, none', add=1)
                    continue
                # note down the regulator type (if given) for the regulator identifier(s) for later processing
                if r_type:
                    for regulator in regulators:
                        self.regulation_type[regulator].add(r_type)

                # TARGET EXTRACTION
                # case: some files contain one column with the regulator identifier and the remaining columns each
                # represent  a target gene identifier => t_col is the start column
                if remaining_cols:
                    targets = {self.get_identifier(source, target) for target in line[t_col:]}
                # case: distinctly specified target column(s)
                else:
                    # process the first target column, which sometimes contains multiple identifiers separated by
                    # forward slashes or commas, and are sometimes surrounded by quotation marks
                    if ' / ' in line[t_col]:
                        targets = {target.strip() for target in line[t_col].replace('"', '').split(' / ')}
                    else:
                        targets = {target.strip() for target in line[t_col].replace('"', '').split(',')}
                    # if more than one target column is given, process those as well and consider them as a possible
                    # mapping to the identifier(s) in the first target column
                    if t_col_secondary:
                        targets = {self.get_identifier(source, target, line[secondary])
                                   for target in targets for secondary in t_col_secondary}
                    # if only one target column is given, process the identifier(s) individually
                    else:
                        targets = {self.get_identifier(source, target) for target in targets}

                # make sure at least one valid target identifier is included in the line
                if not targets:
                    self.log('Regulation initial', 'loading', 'target, none', add=1)

                # SCORE EXTRACTION
                # each regulator -> target interaction has a score between 0.0 and 1.0, where low scores mean that the
                # interaction is tentative and less reliable. most sources do not offer an interaction score and their
                # interactions are assigned a score of 1.0 by default
                score = 1.0
                # if the file contains a score column, extract and process the score
                if score_col:
                    self.log('Regulation initial', 'loading', 'score, yes', add=1)
                    score = float(line[score_col]) / score_divisor

                # STORE INTERACTION
                # since it is not certain at this point which of the regulator or target identifiers are going to be
                # the ones used as a gene symbol in the final mapping, note down an interaction between each regulator
                # and each target identifier and annotate it with the evidence, regulator and target type, source and
                # the score. this is going to be processed, cleaned up and resolved once all regulatory data has
                # been extracted and the identifier mapping is finalised.
                for regulator in regulators:
                    # just to be sure, skip potentially invalid identifiers
                    if not regulator:
                        self.log('Regulation initial', 'loading', 'regulator, empty', add=1)
                        continue
                    # add the regulator gene type
                    if r_type:
                        self.regulation_type[regulator].add(r_type)
                    for target in targets:
                        # just to be sure, skip potentially invalid identifiers
                        if not target:
                            self.log('Regulation initial', 'loading', 'target, empty', add=1)
                            continue
                        # add the target gene type
                        if t_type:
                            self.regulation_type[target].add(t_type)
                        # add the interaction
                        self.regulation_raw[(regulator, target)].sources.add(source)
                        # add the evidence type(s)
                        if evidence == 'predicted':
                            self.regulation_raw[(regulator, target)].predicted = True
                        if evidence == 'experimental':
                            self.regulation_raw[(regulator, target)].experimental = True
                        # add the score
                        if score_col:
                            self.regulation_score[(regulator, target)].add(score)

    def process_regulation_mapping(self):
        """
        Some regulatory data files contain gene aliases for the regulator and/or target identifier. After all regulatory
        data has been extracted from all files, these identifier -> alias mappings contained in these files are
        processed by this function and added to the identifier mapping.
        """
        # initial logging
        self.log('Regulation initial', 'identifier extraction', 'single, to integrate', add=len(self.regulation_single))
        self.log('Regulation initial', 'identifier extraction', 'pairs to integrate', add=len(self.regulation_mapping))

        # go over all identifier pairs and build an index to see how they are interconnected
        by_id = defaultdict(set)                                                        # type: dict[str, set[str]]
        for first, second in self.regulation_mapping.keys():
            by_id[first].add(second)
            by_id[second].add(first)

        # save the current number of identifiers and symbols to later log the changes from adding regulatory mapping
        n_identifiers = len(self.mapping_raw)
        n_symbols = len(self.symbols_raw)
        n_ambiguity = sum(1 for v in self.mapping_raw.values() if len(v) > 1)

        # MAPPING PROCESSING
        for pair, sources in sorted(self.regulation_mapping.items()):
            first, second = pair

            # CASE: only one of them is a known gene symbol from the initial mapping
            # switch around if it's the second one
            if second in self.symbols_main and first not in self.symbols_main:
                first, second = (second, first)

            if first in self.symbols_main and second not in self.symbols_main:
                self.add_raw_identifier(second, first, sources, regulation=True)
                self.log('Regulation inference', 'added pair', 'one is a symbol', add=1)
                continue

            # CASE: only one of them is already in the mapping
            # switch around if it's the second one that is in the mapping
            if second in self.mapping_raw.keys() and first not in self.mapping_raw.keys():
                first, second = (second, first)

            if first in self.mapping_raw.keys() and second not in self.mapping_raw.keys():
                for i in self.mapping_raw[first].keys():
                    self.add_raw_identifier(second, i, sources, regulation=True)
                self.log('Regulation inference', 'added pair', 'one is in initial mapping', add=1)
                continue

            # CASE: both are gene symbols
            if first in self.symbols_main and second in self.symbols_main:
                self.log('Regulation inference', 'added pair', 'both are symbols', add=1)
                continue

            # CASE: both are in the mapping => update the gene symbols they have in common
            if first in self.mapping_raw.keys() and second in self.mapping_raw.keys():
                for i in set(self.mapping_raw[first].keys()).intersection(self.mapping_raw[second].keys()):
                    self.add_raw_identifier(first, i, sources, regulation=True)
                    self.add_raw_identifier(second, i, sources, regulation=True)

                self.log('Regulation inference', 'added pair', 'both are in main mapping', add=1)
                continue

            # CASE: one of them is mapped to by many (=> likely symbol) and the other only maps to one (=>
            # likely identifier)
            if len(by_id[second]) > 1 and len(by_id[first]) == 1:
                first, second = (second, first)

            if len(by_id[first]) > 1 and len(by_id[second]) == 1:
                self.add_raw_identifier(second, first, sources, regulation=True)
                self.log('Regulation inference', 'added pair', 'many and one', add=1)
                continue

            # CASE: both only refer to each other
            if len(by_id[first]) == 1 and len(by_id[second]) == 1:
                # a numeric identifier is likely not a gene symbol and gene symbols tend to be shorter
                # however, there are exceptions in various ID systems
                if first.isdigit():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'Entrez', add=1)
                elif first.startswith('ENSG') and first[4:].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'ENSG', add=1)
                elif first[:3] in ['NM_', 'XM_'] and first[3:].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'NM_, XM_', add=1)
                elif first.startswith('CHEMBL') and first[6:].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'CHEMBL', add=1)
                elif second.startswith('LOC') and second[3:].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'LOC', add=1)
                elif first.upper().startswith('DKFZP') and first[5:].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'DKFZP', add=1)
                elif first[:3] in ['UPI', 'URS'] and len(first) == 13 and first[4:5].isdecimal():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'UPI, URS', add=1)
                elif second.startswith('C') and 'ORF' in second and all(
                        x.isdecimal() for x in second[1:].split('ORF')):
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'ORF', add=1)
                elif len(second) < len(first) and not second.isdigit():
                    first, second = (second, first)
                    self.log('Regulation inference', 'switch (one to one mapping)', 'length', add=1)
                self.add_raw_identifier(second, first, sources, regulation=True)
                self.log('Regulation inference', 'added pair', 'one to one mapping', add=1)
                continue
            self.log('Regulation inference', 'added pair', 'remaining', add=1)

        # log the changes
        self.log('Regulation inference', 'summary', 'identifiers, new from pairs',
                 add=len(self.mapping_raw) - n_identifiers)
        self.log('Regulation inference', 'summary', 'identifiers, new ambiguity from pairs',
                 add=sum(1 for v in self.mapping_raw.values() if len(v) > 1) - n_ambiguity)
        self.log('Regulation inference', 'summary', 'symbols, new from pairs',
                 add=len(self.symbols_raw) - n_symbols)

        # save the current number of identifiers and symbols to later log the changes from adding single identifiers
        n_identifiers = len(self.mapping_raw)
        n_symbols = len(self.symbols_raw)

        # SINGLE IDENTIFIER PROCESSING
        self.log('Regulation inference', 'single mapping', 'total', add=len(self.regulation_single))
        for i, sources in self.regulation_single.items():
            if i in self.mapping_raw.keys():
                self.log('Regulation inference', 'single mapping', 'resolved', add=1)
            elif i in self.mapping_single:
                self.log('Regulation inference', 'single mapping', 'added to single mapping in the meantime', add=1)
            else:
                self.mapping_single[i].update(sources)
                self.log('Regulation inference', 'single mapping', 'added remaining', add=1)

        # log the changes
        self.log('Regulation inference', 'summary', 'identifiers, new from single',
                 add=len(self.mapping_raw) - n_identifiers)
        self.log('Regulation inference', 'summary', 'symbols, new from single',
                 add=len(self.symbols_raw) - n_symbols)

    def add_to_mapping(self, identifier: str, symbol_name: str, sources: set[str]):
        """
        Adds an identifier -> gene symbol mapping.
        :param identifier: name of the identifier
        :param symbol_name: name of the gene symbol that the identifier is mapped to
        :param sources: set of sources supporting the mapping
        """
        # if the identifier is not in the mapping yet, create a new symbol association backed by the given sources
        if identifier not in self.mapping.keys():
            self.mapping[identifier] = Symbol(name=symbol_name, sources=sources)
        symbol = self.mapping[identifier]
        # if the gene symbol has not been added yet, create a new Gene object for it
        if symbol.name not in self.genes.keys():
            self.genes[symbol.name] = Gene()
        # update the gene symbol information
        self.genes[symbol.name].name = symbol.name
        self.genes[symbol.name].aliases[identifier].update(symbol.sources)
        # if the identifier is an Entrez ID, add it to the list of Entrez identifiers of the gene symbol
        if identifier.isdigit():
            self.genes[symbol.name].entrez.add(identifier)

    @staticmethod
    def generate_priority_list(source_labels: list[str], start=4) -> list[tuple[int, str, tuple[bool]]]:
        """
        There are different type of mapping sources, some of which may be more reliable or preferred than others. This
        function generates a list of source combinations in descending order of their reliability.
        :param source_labels: list of source types in order of their reliability, e.g. [A, B, C, D, E, F]
        :param start: where to start assigning the numerical priority
        :return: list of (priority, label, source signature) tuples
        """
        # generate a list of all possible combinations of sources in order of reliability, for example:
        # [(A, B, C, D, E, F), (A, B, C, D, E),..., (B, C, D, E, F),..., (A, B),...,(E, F), (A),..., (F)]
        combs = [tuple(source_labels)]                                                      # type: list[tuple[str]]
        for i in range(len(source_labels) - 1, 0, -1):
            combs += list(combinations(source_labels, i))
        # process the combinations of sources by using the order to assign the priority, compiling the label and
        # the source signature of the combination
        # example: combs = [('A', 'E'), ('A', 'F')] becomes
        # priorities = [(4, 'A + E', (True, False, False, False, True, False),
        #               (5, 'A + F', (True, False, False, False, False, True)]
        priorities = [(i + start, ' + '.join(comb), tuple([label in comb for label in source_labels]))
                      for i, comb in enumerate(combs)]                      # type: list[tuple[int, str, tuple[bool]]]
        return priorities

    def disambiguation(self, priority: int, label: str, pattern: tuple[bool], identifier: str,
                       symbols: dict[str, Symbol], signatures: dict[str, tuple[bool]], ) -> bool:
        """
        The identifier -> symbol mappings are compiled from various types of sources, e.g. manual fixes, the Ensembl
        main mapping, additional mapping from Entrez or additional miRNA mapping from miRBase, as well as from some
        regulation inputs.
        In cases where there are multiple potential mappings, i.e. the sources somewhat disagree,
        there can be combinations of sources that are more reliable than other combinations. For example, if a mapping
        is found in the manual fixes, Ensembl, miRBase AND Entrez, and the other mapping is only found as a side note
        in a regulatory database, this would choose the first mapping and discard the second mapping.
        The function that calls this sub-routine goes over the combinations of sources in order of their reliability.
        This sub-routine then checks how many of the ambiguous mappings for the given identifier match the source
        pattern currently under consideration and applies several disambiguation steps to only those matches.
        :param priority: just used for logging
        :param label: just used for logging
        :param pattern: source signature of the current priority
        :param identifier: the identifier that the function attempts to disambiguate
        :param symbols: dictionary of symbols that the identifier potentially maps to
        :param signatures: dictionary of mapping source signatures for each symbol
        :return: True if the mapping was successfully disambiguated, False otherwise
        """
        label = f'{priority}. priority -- {label}'                                                  # type: str
        # check if there are identifier -> symbol mappings that fit the given source pattern
        matches = list({n for n, signature in signatures.items() if signature == pattern})          # type: list[str]

        # if there are no matches, continue with the next source pattern in the priority list
        if not matches:
            self.log('Untangling', label, 'summary, no match', add=1)
            return False

        # if there is at least one match, try to disambiguate the mapping
        self.log('Untangling', label, 'summary, all processed', add=1)

        # if only exactly one match is found, simply choose that match as the correct mapping
        if len(matches) == 1:
            self.add_to_mapping(identifier, symbols[matches[0]].name, symbols[matches[0]].sources)
            self.log('Untangling', label, 'summary, matches, single', add=1)
            return True

        # if there are multiple matches, apply additional disambiguation criteria to the matched mappings
        self.log('Untangling', label,  'summary, matches, multiple', add=1)

        # LONGEST COMMON PREFIX
        # check if any of the symbols have a common prefix with the identifier and if so, choose the mapping with the
        # longest common prefix. the common prefix should be at least 2 characters long to avoid false positives.
        # example: identifier = ABCD, matches = [ABC, XYZ, KLMN] => choose ABCD -> ABC as the correct mapping
        prefix = LongestOverlap(identifier, matches)                                            # type: LongestOverlap
        # case: single match with the longest common prefix => choose that mapping as the correct one
        if prefix.length > 1 and len(prefix.matches) == 1:
            self.add_to_mapping(identifier, symbols[prefix.matches[0]].name, symbols[prefix.matches[0]].sources)
            self.log('Untangling', label, 'resolved, by longest prefix', add=1)
            return True
        # case: multiple matches that share the longest common prefix
        elif prefix.length > 1 and len(prefix.matches) > 1:
            # check if any of them is equal to the longest common prefix (also check for upper case spelling) and if so,
            # choose that match as the correct mapping
            # example: identifier = ABCD, matches = [ABC, ABCEFG, KLMN] => choose ABCD -> ABC as the correct mapping
            matches = list({x for x in prefix.matches if x == prefix.overlap})
            if not matches:
                matches = list({x for x in prefix.matches if x.upper() == prefix.overlap})
            # case: one of the matches is the longest common prefix => choose that mapping as the correct one
            if len(matches) == 1:
                self.add_to_mapping(identifier, symbols[matches[0]].name, symbols[matches[0]].sources)
                self.log('Untangling', label, 'resolved, by match with longest prefix', add=1)
                return True
            # case: none of the matches is the same as the longest common prefix and therefore there are still multiple
            # candidates for the correct mapping
            # => for the subsequent steps narrow the matches down to the ones that at least share the longest common
            # prefix with the identifier
            # example: identifier = ABCD, matches = [ABCHIJ, ABCEFG, KLMN] => matches = [ABCHIJ, ABCEFG]
            matches = prefix.matches

        # SHORTEST COMMON PREFIX
        # if there are multiple matches that share the longest common prefix with the identifier, check if any of them
        # also has the longest common suffix with the identifier
        # example: identifier = ABCD, matches = [ABC1D, ABCEFG] => choose ABCD -> ABC1D as the correct mapping
        if prefix.length > 1 and len(matches) > 1:
            suffix = LongestOverlap(identifier, matches, suffix=True)                           # type: LongestOverlap

            # case: a single symbol has both the longest common prefix and longest common suffix with the identifier
            # => choose that mapping as the correct one
            if suffix.length > 0 and len(suffix.matches) == 1:
                self.add_to_mapping(identifier, symbols[suffix.matches[0]].name, symbols[suffix.matches[0]].sources)
                self.log('Untangling', label, 'resolved, by longest suffix', add=1)
                return True
            # case: there are still multiple matches that share both the longest common prefix and suffix
            # => narrow down the list of matches to those for subsequent steps
            elif suffix.length > 0 and len(suffix.matches) > 1:
                matches = suffix.matches
                # case: the longest common suffix is a number
                # => choose the match with the shortest numerical suffix
                # example: identifier = ABC1, matches = [ABC.1, ABC.21] => choose ABC1 -> ABC.1 as the correct mapping
                if suffix.overlap.isdigit():
                    # compile the length of the numerical suffix for each match
                    shortest_suffix = sorted({(len(''.join(takewhile(str.isdigit, n[::-1]))), n) for n in matches})
                    # case:
                    # => choose that mapping as the correct one
                    if shortest_suffix[0][0] < shortest_suffix[1][0]:
                        self.add_to_mapping(identifier, symbols[shortest_suffix[0][1]].name,
                                            symbols[shortest_suffix[0][1]].sources)
                        self.log('Untangling', label, 'resolved, by shortest numerical suffix', add=1)
                        return True
                    # case: there is still more than one mapping
                    # => further narrow down the list of matches for subsequent disambiguation steps
                    matches = [n for length, n in shortest_suffix if length == shortest_suffix[0][0]]

        # INITIAL SYMBOLS
        # check if any of the remaining matches can be found in the main list of symbols or the already successfully
        # mapped symbols, which ensures a more consistent and slimmed down usage of symbols
        for collection_name, collection in [('main symbols', self.symbols_main), ('mapped symbols', self.genes)]:
            initial_symbols = list({n for n in matches if n in collection})                     # type: list[str]
            # case: only one of the matches is a primary symbol => choose that mapping as the correct one
            if len(initial_symbols) == 1:
                self.add_to_mapping(identifier, symbols[initial_symbols[0]].name, symbols[initial_symbols[0]].sources)
                self.log('Untangling', label, f'resolved, single known symbols in {collection_name}', add=1)
                return True
            # case: there are still multiple possible matches => narrow down the list further for subsequent steps
            if len(initial_symbols) > 1:
                matches = initial_symbols

        # BASE VERSION
        # if there is still ambiguity that could not be disambiguated with the longest common prefixes and suffixes,
        # check if the identifier seems to be a version and if it is, whether the base version is already known
        if '.' in identifier and identifier.split('.')[-1].isdigit():
            base = '.'.join(identifier.split('.')[:-1])                                         # type: str
            # case: the base version is already in the processed mapping => choose the associated symbol
            if base in self.mapping.keys():
                self.add_to_mapping(identifier, self.mapping[base].name, self.mapping[base].sources)
                self.log('Untangling', label, 'resolved, base version in final mapping', add=1)
                return True
            # case: the base version is among the initial symbol list => choose the base version as the symbol
            if base in self.symbols_main:
                self.add_to_mapping(identifier, base, {'base version'})
                self.log('Untangling', label, 'resolved, base version in initial symbol list', add=1)
                return True
            # case: the base version is in the raw mapping with an unambiguous symbol => choose that symbol
            if base in self.mapping_raw.keys() and len(self.mapping_raw[base]) == 1:
                self.add_to_mapping(identifier, list(self.mapping_raw[base].values())[0].name,
                                    list(self.mapping_raw[base].values())[0].sources)
                self.log('Untangling', label, 'resolved, base version unambiguous in raw mapping', add=1)
                return True

        # NUMBER OF SOURCES
        # use the number of sources as a tiebreaker if possible, choosing the one with the most sources
        n_sources = sorted({(len(symbols[n].sources), n) for n in matches}, reverse=True)
        # if there is a single mapping that is supported by more sources than the others, choose that mapping
        if n_sources[0][0] > n_sources[1][0]:
            self.add_to_mapping(identifier, symbols[n_sources[0][1]].name, symbols[n_sources[0][1]].sources)
            self.log('Untangling', label, 'resolved, by number of sources', add=1)
            return True

        # the mapping could not be disambiguated
        self.log('Untangling', label,  'summary, resolved, no', add=1)
        return False

    def add_self_mapping(self, initial=True):
        """
        All gene symbols should also have an entry in the mapping that maps them to themselves so that they can be
        directly retrieved by TFmiR.
        :param initial: just needed for the log label
        """
        # compile for each symbol in which sources it was mentioned
        self_mapping = defaultdict(set)
        for symbol in self.mapping.values():
            # if the symbol is not yet mapped to itself, store it with the sources that mention it in the index
            if symbol.name not in self.mapping.keys():
                self_mapping[symbol.name].update(symbol.sources)
            # otherwise add the sources to the existing self-mapping
            else:
                self.mapping[symbol.name].sources.update(symbol.sources)
        # go over the index to add an entry for each missing self-mapping
        for symbol, sources in self_mapping.items():
            self.add_to_mapping(symbol, symbol, sources)
        # log the number of added self-mappings depending on where in this function it was done
        if initial:
            self.log('Untangling', 'ambiguity overview', 'self-mapping (initial)', add=len(self_mapping))
        else:
            self.log('Untangling', 'ambiguity overview', 'self-mapping (end)', add=len(self_mapping))

    def resolve_mapping(self):
        """
        The raw mapping of identifiers to gene symbols can contain ambiguous entries where an identifier is mapped to
        multiple gene symbols because the mapping sources are themselves ambiguous or because different mapping sources
        disagree. This function constructs a final, unambiguous mapping by disambiguating ambiguous mappings, making
        sure that each symbol receives itself as an identifier, add single identifiers mapped to themselves, and add
        transformed and intermediate identifier versions.
        """
        # PRE-PROCESSING
        # there are different type of mapping sources, some of which may be more reliable or preferred than others
        source_priority = ['manual fix', 'main fix', 'miRBase', 'main', 'Entrez', 'regulation']     # type: list[str]
        # generate a list of source combinations in descending order of their reliability
        priorities = self.generate_priority_list(source_priority)           # type: list[tuple[int, str, tuple[bool]]]

        # RAW MAPPING DISAMBIGUATION
        # process each raw identifier -> gene symbol(s) mapping and try to resolve ambiguous mapping such that each
        # identifier maps only to a single gene symbol
        for name, ids in sorted(self.mapping_raw.items()):
            # just in case there is an entry without an identifier or without associated symbols
            if not name:
                self.log('Untangling', 'ambiguity overview', 'empty name', add=1)
                continue
            if not ids:
                self.log('Untangling', 'ambiguity overview', 'empty symbols', add=1)
                continue
            # if there is only one associated symbol, there is no ambiguity and the mapping can be stored immediately
            if len(ids) == 1:
                self.add_to_mapping(name, list(ids.values())[0].name, list(ids.values())[0].sources)
                self.log('Untangling', 'ambiguity overview', 'unambiguous', add=1)
                continue

            # there are multiple symbols associated with the identifier that needs to be disambiguated
            self.log('Untangling', 'ambiguity overview', 'ambiguous', add=1)

            # priority 1: Entrez identifiers should not be gene symbols if possible
            if name.isdigit():
                self.log('Untangling', '1. priority -- Entrez', 'all', add=1)
                # get the non-Entrez symbols mapped to this identifier
                symbols_without_entrez = {s for s in ids.keys() if
                                          not self.symbols_raw[s].get_Entrez_ids(ignore={name})}
                # if there is only one non-Entrez symbol, map the identifier to that symbol
                if len(symbols_without_entrez) == 1:
                    self.add_to_mapping(name, ids[list(symbols_without_entrez)[0]].name,
                                        ids[list(symbols_without_entrez)[0]].sources)
                    self.log('Untangling', '1. priority -- Entrez', 'single non-numerical symbol', add=1)
                    continue
                # otherwise simply log additional information
                elif not symbols_without_entrez:
                    self.log('Untangling', '1. priority -- Entrez', 'no non-numerical symbol', add=1)
                else:
                    self.log('Untangling', '1. priority -- Entrez', 'multiple non-numerical symbols', add=1)
            else:
                self.log('Untangling', '1. priority -- Entrez', 'not Entrez ID', add=1)

            # priority 2: mappings from the manual fixes file
            manual_fix = list({name for name, info in ids.items() if info.in_manual_fix})
            # if there is a unique match, choose that symbol as the correct one
            if len(manual_fix) == 1:
                self.add_to_mapping(name, ids[manual_fix[0]].name, ids[manual_fix[0]].sources)
                self.log('Untangling', '2. priority -- manual fix', 'single match', add=1)
                continue
            # otherwise simply log additional information
            elif len(manual_fix) > 1:
                self.log('Untangling', '2. priority -- manual fix', 'multiple', add=1)
            elif not manual_fix:
                self.log('Untangling', '2. priority -- manual fix', 'no match', add=1)

            # priority 3: choose the identifier itself if it is amongst the symbols
            if name in ids:
                self.add_to_mapping(name, ids[name].name, ids[name].sources)
                self.log('Untangling', '3. priority -- identifier amongst symbols', 'yes', add=1)
                continue
            self.log('Untangling', '3. priority -- identifier amongst symbols', 'no', add=1)

            # priorities 4+: prioritise symbols based on the type of source
            # if none of the previous steps worked or were applicable, go over the source priority list
            # starting with the most reliable or preferred combination of sources and try to disambiguate the mapping
            # based only on the symbols matching the source pattern
            success = False
            # compile the source type pattern for each symbol associated with the current identifier
            # (the order of this should match the source_priority list at the beginning of the function)
            signatures = {n: (s.in_manual_fix, s.in_main_cleanup, s.in_miRBase, s.in_main, s.in_Entrez, s.in_regulation)
                          for n, s in ids.items()}
            # apply the disambiguation steps for each source combination
            for priority, label, pattern in priorities:
                # if successful, the disambiguation function already added the correct mapping and we can stop
                if success := self.disambiguation(priority, label, pattern, name, ids, signatures):
                    break
                # if unsuccessful, move onto the next source combination and try again
            # if none of the disambiguation attempts worked, store the mapping as ambiguous
            if not success:
                self.ambiguous[name] = ids

        # post-process the logged information of the source type-based disambiguation
        for priority, label, pattern in priorities:
            label = f'{priority}. priority -- {label}'
            log = self.get('Untangling', label)
            # remove priority cases from the log that did not occur in any of the ambiguous mappings
            if set(log.keys()) == {'summary, no match'}:
                del self.log_dict['Untangling'][label]
                continue
            # remove the skipped entry from the remaining log entries (originally added to display them in order)
            del self.log_dict['Untangling'][label]['summary, no match']
            # compute for each priority, how many identifiers could be resolved in total
            resolved = self.get('Untangling', label, 'summary, all processed', default=0)
            resolved -= self.get('Untangling', label, 'summary, resolved, no', default=0)
            self.log('Untangling', label, 'summary, resolved, yes', add=resolved)
            #
            self.log('Untangling', 'disambiguation summary', 'z. single signature match',
                     add=self.get('Untangling', label, 'summary, matches, single', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. base version in final mapping',
                     add=self.get('Untangling', label, 'resolved, base version in final mapping', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. by longest prefix',
                     add=self.get('Untangling', label, 'resolved, by longest prefix', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. by longest suffix',
                     add=self.get('Untangling', label, 'resolved, by longest suffix', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. by match with longest prefix',
                     add=self.get('Untangling', label, 'resolved, by match with longest prefix', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. by number of sources',
                     add=self.get('Untangling', label, 'resolved, by number of sources', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. by shortest numerical suffix',
                     add=self.get('Untangling', label, 'resolved, by shortest numerical suffix', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. single known symbols in main symbols',
                     add=self.get('Untangling', label, 'resolved, single known symbols in main symbols', default=0))
            self.log('Untangling', 'disambiguation summary', 'z. single known symbols in mapped symbols',
                     add=self.get('Untangling', label, 'resolved, single known symbols in mapped symbols', default=0))

        # compile and log some disambiguation statistics
        self.log('Untangling', 'disambiguation summary', 'all examined', 
                 add=self.get('Untangling', 'ambiguity overview', 'ambiguous', default=0))
        self.log('Untangling', 'disambiguation summary', 'resolved, no', add=len(self.ambiguous))
        self.log('Untangling', 'disambiguation summary', 'resolved, yes',
                 add=self.get('Untangling', 'disambiguation summary', 'all examined') - len(self.ambiguous))

        self.log('Untangling', 'disambiguation summary', 'z. Entrez: single non-numerical symbol',
                 add=self.get('Untangling', '1. priority -- Entrez', 'single non-numerical symbol', default=0))
        self.log('Untangling', 'disambiguation summary', 'z. manual fix',
                 add=self.get('Untangling', '2. priority -- manual fix', 'single match', default=0))
        self.log('Untangling', 'disambiguation summary', 'z. identifier amongst symbols',
                 add=self.get('Untangling', '3. priority -- identifier amongst symbols', 'yes', default=0))

        # SELF MAPPING
        # all gene symbols should also have an entry in the mapping that maps them to themselves so that they can be
        # directly retrieved by TFmiR
        self.add_self_mapping()

        # AMBIGUOUS IDs
        # go over the not yet disambiguated identifier -> gene symbol mappings and map the identifiers and symbols to
        # themselves so that they are not lost but also no potentially false mappings are added
        for identifier, symbols in self.ambiguous.items():
            self.log('Untangling', 'self-mapping (ambiguous)', 'ambiguous identifiers', add=1)
            # merge the sources from all mappings involving the identifier als sources for the self-mapping
            sources = set()
            for name, symbol in symbols.items():
                sources.update(symbol.sources)
            # if the identifier is not yet included in the mapping from previous steps, add it as a self-mapping
            if identifier not in self.mapping.keys():
                self.add_to_mapping(identifier, identifier, sources)
                self.log('Untangling', 'self-mapping (ambiguous)', 'ambiguous identifier added to mapping', add=1)
            else:
                self.log('Untangling', 'self-mapping (ambiguous)', 'identifier already included in mapping', add=1)

        # SINGLE IDs
        # finally, add all identifiers that were at various points read as input without an accompanying gene symbol
        for identifier, sources in self.mapping_single.items():
            self.log('Untangling', 'self-mapping (single)', 'total processed', add=1)
            # if any of the other information or previous steps already mapped the identifier to a gene symbol,
            # there is nothing more do to for the identifier
            if identifier in self.mapping.keys():
                self.log('Untangling', 'self-mapping (single)', 'already in mapping', add=1)
                continue
            # otherwise map the identifier to itself with all the sources that mentioned the identifier
            self.log('Untangling', 'self-mapping (single)', 'not yet in mapping', add=1)
            self.add_to_mapping(identifier, identifier, sources)

        # TRANSFORMED IDs
        # some identifiers were processed and transformed to fix spelling mistakes or achieve a consistent spelling
        # of the same identifier. this step adds these transformed identifiers to the mapping as well
        self.log('Untangling', 'transformed ID processing', 'processed identifier, all',
                 add=len(self.transformed))
        # first, turn the final_identifier: [raw identifiers,...]into a raw_identifier: [final_identifiers,
        # ...] index
        transformed_index = defaultdict(set)
        for identifier, raw in self.transformed.items():
            # log whether the final identifier is already included in the mapping
            if identifier in self.mapping.keys():
                self.log('Untangling', 'transformed ID processing', 'processed identifier, already in mapping',
                         add=1)
            else:
                self.log('Untangling', 'transformed ID processing', 'processed identifier, not yet in mapping',
                         add=1)
            # flip the index
            for t in raw:
                transformed_index[t].add(identifier)
        # log how many unique raw identifiers have to be processed
        self.log('Untangling', 'transformed ID processing', 'raw identifiers, total', add=len(transformed_index))

        # try to add each raw identifier to the mapping
        for raw, identifiers in transformed_index.items():
            # just in case there are somehow raw identifiers without associated final identifiers
            if not identifiers:
                self.log('Untangling', 'transformed ID processing', 'raw identifiers, no identifier', add=1)
                continue
            # go over the associated final identifiers and retrieve the corresponding gene symbols from the mapping
            symbols = list({self.mapping[i].name for i in identifiers if i in self.mapping.keys()})
            # log whether the raw identifier is uniquely or ambiguously associated with final identifiers and
            # whether
            # the corresponding symbols are unique or ambiguous
            if len(identifiers) > 1:
                self.log('Untangling', 'transformed ID processing',
                         'raw identifiers, more than one final identifier', add=1)
                if not symbols:
                    self.log('Untangling', 'transformed ID processing',
                             'raw identifiers, more than one final identifier - none known', add=1)
                elif len(symbols) > 1:
                    self.log('Untangling', 'transformed ID processing',
                             'raw identifiers, more than one final identifier - ambiguous symbols', add=1)
                    continue
                self.log('Untangling', 'transformed ID processing',
                         'raw identifiers, more than one final identifier - unambiguous symbols', add=1)
            if len(identifiers) == 1:
                self.log('Untangling', 'transformed ID processing', 'raw identifiers, single final identifier',
                         add=1)

                if not symbols:
                    self.log('Untangling', 'transformed ID processing',
                             'raw identifiers, single final identifier - not known', add=1)
                else:
                    self.log('Untangling', 'transformed ID processing',
                             'raw identifiers, single final identifier - known', add=1)
            # if neither the final identifier(s) nor the raw identifier are known, map the raw identifier to itself
            if not symbols and raw not in self.mapping.keys():
                self.add_to_mapping(raw, raw, set())
                self.log('Untangling', 'transformed ID processing', 'mapping, all unknown', add=1)
                continue
            # if the final identifier(s) are not in the mapping but the raw identifier is, map it and the final
            # identifier(s) to the associated symbol
            elif not symbols and raw in self.mapping.keys():
                self.log('Untangling', 'transformed ID processing', 'mapping, only raw identifier known', add=1)
                symbols = [self.mapping[raw].name]
            # if the final identifier(s) are in the mapping, simply proceed with mapping the raw to symbol
            # associated
            # with the final identifier(s)
            elif symbols and raw not in self.mapping.keys():
                self.log('Untangling', 'transformed ID processing', 'mapping, only final identifier known', add=1)
            # if both the final identifier(s) and the raw identifier are in the mapping, make sure they map to
            # the same
            # gene symbol and if not, log the incident
            else:
                self.log('Untangling', 'transformed ID processing', 'mapping, both known', add=1)
                if self.mapping[raw].name not in symbols:
                    self.log('Untangling', 'transformed ID processing', 'mapping, both known - ambiguous', add=1)
                    continue
                self.log('Untangling', 'transformed ID processing', 'mapping, both known - unambiguous', add=1)

            self.log('Untangling', 'transformed ID processing', 'unambiguous mapping, all', add=1)
            # the previous steps already made sure that the symbol list only has one element
            symbol = symbols[0]
            # add the raw identifier to the identifier set, which are then all mapped to the symbol
            identifiers.add(raw)
            # remove all identifiers from the set that are already included in the mapping
            identifiers = {i for i in identifiers if
                           i and i != symbol and i not in self.mapping.keys() and i not in self.ambiguous.items()}
            # log the case where all identifiers were already included in the mapping
            if not identifiers:
                self.log('Untangling', 'transformed ID processing', 'unambiguous mapping, all already known', add=1)
                continue
            self.log('Untangling', 'transformed ID processing', 'unambiguous mapping, some unknown', add=1)
            # map all identifiers to the symbol
            for identifier in identifiers:
                self.add_to_mapping(identifier, symbol, {'in-house processing'})
                self.log('Untangling', 'transformed ID processing', 'unambiguous mapping, new identifiers', add=1)

        # SELF-MAPPING
        # run the self-mapping again just to make sure that all symbols that may have been added in the meantime are
        # also represented in the mapping
        self.add_self_mapping()

    def resolve_regulation(self):
        """
        Now that the final identifier to gene symbol mapping is completed, go over the collected interactions and
        store them by mapping the regulator and target identifiers to the corresponding symbols, which reduces
        interaction duplicates and unifies the available data for the interactions. Then add scores to the interactions.
        """
        # INTERACTIONS
        # use the identifier -> gene symbol mapping to process and store regulator-target interactions
        for (regulator, target), info in self.regulation_raw.items():
            # make sure the regulator and target have known identifiers and if not, log the problem
            if regulator not in self.mapping.keys() or target not in self.mapping.keys():
                self.log('Resolve regulation', 'interactions, skipped', 'total', add=1)
                for label, identifier in {'regulator': regulator, 'target': target}.items():
                    if identifier in self.mapping.keys():
                        continue
                    self.log('Resolve regulation', 'interactions, skipped', f'missing {label}', add=1)
                    self.log('Resolve regulation', f'{label}, missing', 'total', update=identifier)
                    if identifier in self.mapping_raw.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in raw mapping', update=identifier)
                    elif identifier in self.ambiguous.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in ambiguous', update=identifier)
                    elif identifier in self.mapping_single.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in single', update=identifier)
                    elif identifier in self.transformed.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in transformed', update=identifier)
                continue
            # map the regulator and target identifier to the corresponding gene symbols
            regulator = self.mapping[regulator].name            # type: str
            target = self.mapping[target].name                  # type: str
            # store and log some additional information about the interaction
            if regulator == target:
                self.log('Resolve regulation', 'interactions, added', 'self-regulation', add=1)
            if info.experimental:
                self.regulation[(regulator, target)].experimental = True
                self.log('Resolve regulation', 'interactions, added', 'experimental', add=1)
            if info.predicted:
                self.regulation[(regulator, target)].predicted = True
                self.log('Resolve regulation', 'interactions, added', 'predicted', add=1)
            if info.predicted and info.experimental:
                self.log('Resolve regulation', 'interactions, added', 'both', add=1)
            # update the sources of the interaction and set the interaction type based on the regulator and target type
            self.regulation[(regulator, target)].sources.update(info.sources)
            self.regulation[(regulator, target)].type = f'{self.genes[regulator].type}-{self.genes[target].type}'
            # note down in the regulator and target genes that they are involved in known interactions
            self.genes[regulator].in_regulation = True
            self.genes[regulator].is_regulator = True
            self.genes[target].in_regulation = True
            self.genes[target].is_target = True

            self.log('Resolve regulation', 'interactions, added', 'total', add=1)

        # SCORES
        # map the identifier based interaction score information to symbols
        regulation_scores = defaultdict(set)                                # type: dict[tuple[str, str], set[float]]
        for (regulator, target), scores in self.regulation_score.items():
            # make sure the regulator and target have known identifiers and if not, log the problem
            if regulator not in self.mapping.keys() or target not in self.mapping.keys():
                self.log('Resolve regulation', 'scores, skipped', 'total', add=1)
                for label, identifier in {'regulator': regulator, 'target': target}.items():
                    if identifier in self.mapping.keys():
                        continue
                    self.log('Resolve regulation', 'scores, skipped', f'missing {label}', add=1)
                    self.log('Resolve regulation', f'{label}, missing', 'total', update=identifier)
                    if identifier in self.mapping_raw.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in raw mapping', update=identifier)
                    if identifier in self.ambiguous.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in ambiguous', update=identifier)
                    if identifier in self.mapping_single.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in single', update=identifier)
                    if identifier in self.transformed.keys():
                        self.log('Resolve regulation', f'{label}, missing', 'in transformed', update=identifier)
                continue
            # map the regulator and target identifiers to the corresponding symbols
            regulator = self.mapping[regulator].name
            target = self.mapping[target].name
            # add the score to the interaction
            regulation_scores[(regulator, target)].update(scores)

        # add scores to interactions for which they were provided
        for (regulator, target), scores in regulation_scores.items():
            # make sure that the interaction is known and that there is score information given
            if (regulator, target) not in self.regulation.keys():
                self.log('Resolve regulation', 'scores, skipped', 'unknown regulation', add=1)
                continue
            if not scores:
                self.log('Resolve regulation', 'scores, skipped', 'no scores', add=1)
                continue
            self.log('Resolve regulation', 'scores, added', 'all', add=1)
            # if there is conflicting scores given for the interaction, choose the highest one
            if len(scores) > 1:
                self.log('Resolve regulation', 'scores, added', 'ambiguous', add=1)
                score = max(scores)
            else:
                self.log('Resolve regulation', 'scores, added', 'unambiguous', add=1)
                score = list(scores)[0]
            # log some statistics about the score
            if score == 1.0:
                self.log('Resolve regulation', 'scores, added', 'value, 1.0', add=1)
            elif score > 1.0:
                self.log('Resolve regulation', 'scores, added', 'value, > 1.0', add=1)
            elif score < 0.0:
                self.log('Resolve regulation', 'scores, added', 'value, < 0', add=1)
            else:
                self.log('Resolve regulation', 'scores, added', 'value, between 0 and 1.0', add=1)
            # add the score to the corresponding regulation and make sure that it is between 0 and 1
            self.regulation[(regulator, target)].score = max(min(score, 1.0), 0.0)

        # PPI only
        print(len(self.regulation))
        print(self.ppi_only)
        for regulation in self.regulation.values():
            if len(self.ppi_only.intersection(regulation.sources)) == len(regulation.sources):
                regulation.ppi_only = True

    # ANNOTATIONS ------------------------------------------------------------------------------------------------------
    def load_annotations(self):
        """
        Annotate genes with gene type, tissue, process and disease information. This should be done after the final
        mapping has been generated and fully processed.
        """
        # identify TFs and miRNAs
        if self.exec('load gene types', self.load_gene_types, flag=self.cfg.run.skip_gene_type, is_sub=True):
            return True
        # annotate miRNAs with disease spectrum width information
        if self.exec('load DSW', self.load_DSW, flag=self.cfg.run.skip_DSW, is_sub=True):
            return True

        # annotate genes with disease, process and tissue information
        def extract():
            for a in self.cfg.paths.annotations:
                self.load_annotation(file_path=a.path, annotation_type=a.type, gene_col=a.gene,
                                     annotation_col=a.annotation, delimiter=a.delimiter, header=a.header,
                                     gene_remaining_cols=a.remaining_cols, split_annotation=a.split,
                                     min_cols=a.min_cols)
        if self.exec('extract annotations', extract, flag=self.cfg.run.skip_annotations, is_sub=True):
            return True

        # note if no annotations could be extracted even though annotation loading was enabled
        if not self.cfg.run.skip_annotations and not self.annotations_raw:
            self.cfg.run.skip_annotations = 'no annotations extracted'
        # process and merge the raw annotations
        if self.exec('process annotations', self.annotation_merge, flag=self.cfg.run.skip_annotations, is_sub=True):
            return True

        return False

    def load_gene_types(self):
        """
        Load gene type information (TF, miRNA) from the specified file and integrate it with type information extracted
        from regulatory data.
        """
        type_dict = defaultdict(set)                                                    # type: dict[str, set[str]]

        # FROM FILE
        # load gene type information from the specified file
        load_dict = defaultdict(set)                                                    # type: dict[str, set[str]]
        with open(self.cfg.paths.gene_type.path, 'r') as file:
            for line in csv.DictReader(file, delimiter='\t'):
                identifier, _ = self.process_identifier(line['gene'])
                # skip unknown identifiers
                if identifier not in self.mapping:
                    self.log('Annotations', 'gene types, load', 'identifiers, not known', update=identifier)
                    continue
                # skip lines without gene type information
                gene_type = line['gene type'].strip()                                   # type: str
                if not gene_type:
                    self.log('Annotations', 'gene types, load', 'types, empty', add=1)
                    continue
                # skip lines with unknown gene types
                if gene_type not in ['TF', 'miRNA', 'gene']:
                    self.log('Annotations', 'gene types, load', 'types, unknown', add=1)
                    continue
                self.log('Annotations', 'gene types, load', 'identifiers, known', update=identifier)
                load_dict[self.mapping[identifier].name].add(gene_type)
        self.log('Annotations', 'gene types, load', 'identifiers, recognized symbols', add=len(load_dict))
        # add the type information
        for symbol, types in load_dict.items():
            if len(types) > 1:
                self.log('Annotations', 'gene types, ambiguous', 'all, from loading', update=symbol)
                self.log('Annotations', 'gene types, load', 'types, ambiguous', update=symbol)
            else:
                self.log('Annotations', 'gene types, load', f'types, {list(types)[0]}', update=symbol)
            type_dict[symbol].update(types)

        # FROM REGULATION
        # genes already annotated with type information from the file
        already_annotated = set(type_dict.keys())                                           # type: set[str]
        # map the identifiers annotated with a gene type by regulatory information to symbols
        regulation_types = defaultdict(set)                                                 # type: dict[str, set[str]]
        for identifier, types in self.regulation_type.items():
            # skip unknown identifiers, though this should not happen
            if identifier not in self.mapping:
                self.log('Annotations', 'gene types, regulation', 'identifier, unknown', update=identifier)
                continue
            regulation_types[self.mapping[identifier].name].update(types)
        # process type information obtained during regulation processing
        for symbol, types in regulation_types.items():
            # log if the type information is inconsistent
            if len(types) > 1:
                self.log('Annotations', 'gene types, ambiguous', 'all, from regulation', add=1)
                self.log('Annotations', 'gene types, regulation', 'types, ambiguous', add=1)
            else:
                self.log('Annotations', 'gene types, regulation', f'types, {list(types)[0]}', update=symbol)
            # log if the type annotation from regulation is consistent with that from the file
            if symbol in already_annotated:
                self.log('Annotations', 'gene types, regulation', 'already annotated, all', update=symbol)
                if types == type_dict[symbol]:
                    self.log('Annotations', 'gene types, regulation', 'already annotated, match', update=symbol)
                else:
                    self.log('Annotations', 'gene types, regulation', 'already annotated, mismatch', update=symbol)
            else:
                self.log('Annotations', 'gene types, regulation', 'already annotated, no', update=symbol)
            # add the type information
            type_dict[symbol].update(types)

        # process the combined gene type information and add it to the respective genes if it is consistent
        for symbol, types in type_dict.items():
            # try to figure out the correct type if there is inconsistent type information
            if len(types) > 1:
                self.log('Annotations', 'gene types, ambiguous', 'all', add=1)
                # use the identifier type to see if the gene is a miRNA
                if 'miRNA' in types and is_mirna_symbol(symbol):
                    self.log('Annotations', 'gene types, ambiguous', 'fixed, miRNA symbol', add=1)
                    types = {'miRNA'}
                # if the type information is still inconsistent, do not add anything
                else:
                    self.log('Annotations', 'gene types, ambiguous', 'unresolved', add=1)
                    continue
            # in some data sets miRNA host genes are mistakenly labelled as miRNAs
            if types == {'miRNA'} and symbol.lower().startswith('mir') and symbol.lower().endswith('hg'):
                self.log('Annotations', 'gene types, processing', 'miRNA host gene', add=1)
                continue
            # if the type information is consistent or could be made consistent, add it to the gene
            self.genes[symbol].type = list(types)[0]
            self.log('Annotations', 'gene types, final summary', f'{list(types)[0]}s', update=symbol)

        # see if there are genes with miRNA IDs or accessions that are not labelled as miRNAs
        for symbol, gene in self.genes.items():
            # skip genes that have already been annotated
            if gene.type != 'gene':
                continue
            self.log('Annotations', 'gene types, inference', 'all', add=1)
            if is_mirna_symbol(symbol) and not symbol.endswith('HG'):
                self.log('Annotations', 'gene types, inference', 'miRNA symbol', add=1)
                gene.type = 'miRNA'
            elif is_miRBase_ID(symbol, self.cfg.miRNA_prefix):
                self.log('Annotations', 'gene types, inference', 'miRBase ID', add=1)
            elif is_miRBase_accession(symbol):
                self.log('Annotations', 'gene types, inference', 'miRBase accession', add=1)

        # also log the number of normal genes
        self.log('Annotations', 'gene types, final summary', 'genes',
                 add=sum(1 for gene in self.genes.values() if gene.type == 'gene'))

    def load_DSW(self):
        """
        Load information about the disease spectrum width (DSW) of miRNAs.
        """
        with open(self.cfg.paths.DSW.path, 'r') as file:
            for line in file:
                cols = line.strip().split('\t')                             # type: list[str]
                # skip rows with missing columns and the header row
                if len(cols) != 3 or cols == ['miRNA', 'n', 'DSW']:
                    continue
                miRNA, _, DSW = cols                                        
                # use the mapping to find and annotate the corresponding miRNA gene
                if miRNA in self.mapping:
                    self.genes[self.mapping[miRNA].name].dsw = float(DSW)
                    self.log('Annotations', 'DSW', 'annotated miRNAs', update=self.mapping[miRNA].name)
                else:
                    self.log('Annotations', 'DSW', 'unknown identifier', update=miRNA)

    def process_annotation(self, annotation: str, annotation_type: str):
        """
        Processes an annotation to remove unnecessary information and ensure consistent spelling.
        :param annotation: annotation string to process
        :param annotation_type: should be either tissue, disease or process
        :return: fixed annotation
        """
        input_annotation = annotation           # type: str

        # INITIAL ADJUSTMENTS
        # remove excessive whitespace and cast to lower case for consistent spelling
        annotation = ' '.join(annotation.lower().strip().split())
        # remove unnecessary information
        if '(' in annotation and ')' in annotation:
            # Keratinocyte Proliferation(26402295;26402295;25654102)
            # Immune System(Xiao's Cell2010)
            # Anti-Cell Proliferation(Hwang Etal Bjc2007)
            for stripped in annotation.split('(')[1].strip(')').split(';'):
                if stripped.isdigit() or stripped.startswith('xiao') or stripped.startswith('hwang'):
                    annotation = annotation.split('(')[0].strip()
                    break

        # REMOVAL and REPLACEMENT
        # remove undesirable components
        for to_remove in self.cfg.format.remove:
            annotation = annotation.replace(to_remove, '')
        # remove surrounding whitespace generated by some fixes
        annotation = ' '.join(annotation.strip().split())
        # replace larger segments
        for to_replace, replacement in self.cfg.format.replace_multiple.items():
            annotation = annotation.replace(to_replace, replacement)
        # remove surrounding whitespace generated by some fixes
        annotation = ' '.join(annotation.strip().split())
        # adjust individual words in the annotation
        for replacement, pattern in self.cfg.format.replace_single.items():
            annotation = pattern.sub(replacement, annotation)
        # remove surrounding whitespace generated by some fixes
        annotation = ' '.join(annotation.strip().split())
        for replacement, pattern in self.cfg.format.replace_regex.items():
            annotation = pattern.sub(replacement, annotation)
        # remove surrounding whitespace generated by some fixes
        annotation = ' '.join(annotation.strip().split())

        # CASE ADJUSTMENT
        for pattern in self.cfg.format.upper_complex:
            annotation = re.sub(pattern, to_upper, annotation)
        for replacement, pattern in self.cfg.format.upper_simple.items():
            annotation = pattern.sub(replacement, annotation)
        for replacement, pattern in self.cfg.format.capitalize.items():
            annotation = pattern.sub(replacement, annotation)
        for replacement, pattern in self.cfg.format.mixed.items():
            annotation = pattern.sub(replacement, annotation)
        for replacement, pattern in self.cfg.format.lower.items():
            annotation = pattern.sub(replacement, annotation)

        # check if the annotation is invalid
        if not annotation or annotation.lower() in self.cfg.format.invalid:
            self.log('Annotations', 'annotations, fix', 'all, invalid', add=1)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, invalid', add=1)
            return ''

        # LOGGING
        if annotation != input_annotation:
            self.log('Annotations', 'annotations, fix', 'all, processed', add=1)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, processed', add=1)
        else:
            self.log('Annotations', 'annotations, fix', 'all, not altered', add=1)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, not altered', add=1)

        return annotation

    def load_annotation(self, file_path: str, annotation_type: str, gene_col: str | int,
                        annotation_col: str | int, header=True, delimiter='\t', gene_remaining_cols=False,
                        min_cols=0, split_annotation=''):
        """
        Processes an annotation file and adds the information to the specified genes.
        :param file_path: path to the annotation file
        :param annotation_type: should be tissue, disease or process
        :param gene_col: name (or index) of the column containing the gene identifier
        :param annotation_col: name (or index) of the column containing annotation(s)
        :param delimiter: string that separates the columns in a row
        :param header: True if the file has column names specified in the first row, False otherwise
        :param gene_remaining_cols: True if the columns after the gene column also contain gene identifiers
        :param min_cols: minimum number of columns a row needs to have to be considered valid for processing
        :param split_annotation: separator by which to split the annotation column, if it contains multiple annotations
        """
        with open(file_path, 'r', encoding='latin-1') as file:
            # if the file has a header, use the column names, otherwise simply split the lines
            reader = csv.DictReader(file, delimiter=delimiter) if header else file

            for line in reader:
                # if the file has a header, use the column names, otherwise simply split the lines and use indices
                cols = line if header else line.strip().split(delimiter)            # type: list[str] | dict[str, str]
                self.log('Annotations', 'annotations, load', 'rows, all', add=1)
                if len(cols) < min_cols:
                    self.log('Annotations', 'annotations, load', 'rows, missing columns', add=1)
                    continue
                if not cols[annotation_col] or not cols[gene_col]:
                    self.log('Annotations', 'annotations, load', 'rows, missing symbol or annotation', add=1)
                    continue
                self.log('Annotations', 'annotations, load', 'rows, full', add=1)
                # extract all identifiers
                if gene_remaining_cols:
                    identifiers = {self.process_identifier(g)[0] for g in cols[gene_col:]}      # type: set[str]
                else:
                    identifiers = {self.process_identifier(cols[gene_col])[0]}
                # extract all annotations and fix them if necessary
                if split_annotation:
                    annotations = {a for a in cols[annotation_col].split(', ')}                 # type: set[str]
                else:
                    annotations = {cols[annotation_col]}
                # make sure at least one identifier and one annotation could be extracted
                if not identifiers:
                    self.log('Annotations', 'annotations, load', 'rows, invalid identifiers', add=1)
                if not annotations:
                    self.log('Annotations', 'annotations, load', 'rows, empty annotations', add=1)
                if not annotations or not identifiers:
                    continue
                # add the annotation(s)
                for identifier in identifiers:
                    self.log('Annotations', 'annotations, load', 'identifiers, all', update=identifier)
                    # the identifier is not in the identifier mapping
                    if identifier not in self.mapping:
                        self.log('Annotations', 'annotations, load', 'identifiers, not known', update=identifier)
                        continue
                    # use the mapping to find the gene symbol associated with the identifier
                    symbol = self.mapping[identifier].name                                              # type: str
                    # log some information
                    self.log('Annotations', 'annotations, load', 'identifiers, known', update=identifier)
                    self.log('Annotations', 'annotations, load', 'identifiers, recognized symbols', update=symbol)
                    # add the symbol to each annotation
                    for annotation in annotations:
                        self.annotations_raw[annotation_type][annotation].add(symbol)                    

    def annotation_merge(self):
        """
        Processes all raw annotations to adjust spelling and the like and merges the information.
        """
        # {annotation type: {symbol,...}}
        symbols_by_type = defaultdict(set)                                              # type: dict[str, set[str]]
        # process and merge the annotation information
        for annotation_type, annotation_dict in self.annotations_raw.items():
            for annotation, symbols in annotation_dict.items():
                # process the annotation
                annotation = self.process_annotation(annotation, annotation_type)
                # invalid annotation
                if not annotation:
                    continue
                symbols_by_type[annotation_type].update(symbols)
                # map the annotation to genes
                self.annotations[annotation_type][annotation].update(symbols)
                for symbol in symbols:
                    if annotation_type == 'process':
                        self.genes[symbol].processes.add(annotation)
                    elif annotation_type == 'disease':
                        self.genes[symbol].diseases.add(annotation)
                    elif annotation_type == 'tissue':
                        self.genes[symbol].tissues.add(annotation)

        # LOGGING
        # log how many symbols have been annotated
        for annotation_type, symbols in symbols_by_type.items():
            self.log('Annotations', f'annotations, {annotation_type}', 'symbols', add=len(symbols))
            
        for annotation_type, raw_dict in self.annotations_raw.items():
            raw = len(raw_dict)
            final = len(self.annotations[annotation_type])
            # log how many unique annotations were extracted from the input
            self.log('Annotations', 'annotations, fix', 'all, input', add=raw)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, input', add=raw)
            # log how many unique annotations exist after processing
            self.log('Annotations', 'annotations, fix', 'all, total', add=final)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, total', add=final)
            self.log('Annotations', f'annotations, {annotation_type}', 'all', add=final)
            # log how many annotations were different spellings of the same term and could be merged
            self.log('Annotations', 'annotations, fix', 'all, reduction', add=raw - final)
            self.log('Annotations', 'annotations, fix', f'{annotation_type}, reduction', add=raw - final)

    # POST-PROCESSING --------------------------------------------------------------------------------------------------
    def final_cleanup(self):
        """
        Move symbols (genes) for which neither annotation nor regulation data exist into a separate table, which
        removes putative and other little supported genes introduced by some sources. This way the information is still
        available if needed without cluttering data retrieval and processing by TFmiR.
        Identifiers in the mapping that map to discarded gene symbols are similarly moved to a separate table to speed
        up identifier mapping and to ensure that all gene symbols in the mapping can be found in the gene table.
        """
        # remove symbols (genes) for which neither annotation nor regulation data exist
        self.log('Final cleanup', 'symbols', 'all', add=len(self.genes))
        for name in list(self.genes.keys()):
            gene = self.genes[name]                 # type: Gene
            if not gene.in_regulation and not gene.tissues and not gene.diseases and \
                    not gene.processes and not gene.dsw:
                self.log('Final cleanup', 'symbols', 'removed', add=1)
                self.genes_discarded[name] = gene
                del self.genes[name]
            else:
                self.log('Final cleanup', 'symbols', 'kept', add=1)

        # remove identifiers from the mapping that map to discarded symbols (genes)
        self.log('Final cleanup', 'identifiers', 'all', add=len(self.mapping))
        for identifier in list(self.mapping.keys()):
            symbol = self.mapping[identifier]       # type: Symbol
            if symbol.name not in self.genes.keys():
                self.log('Final cleanup', 'identifiers', 'removed', add=1)
                self.mapping_discarded[identifier] = symbol
                del self.mapping[identifier]
            else:
                self.log('Final cleanup', 'identifiers', 'kept', add=1)

    # DATABASE GENERATION ----------------------------------------------------------------------------------------------
    def build_database(self) -> bool:
        """
        Generate all database tables and other files that are enabled in the configuration.
        :return: True if an error occurred in a database building step, False otherwise
        """
        if self.exec('generate gene table', self.gene_table, flag=self.cfg.run.skip_gene_table, is_sub=True):
            return True

        if self.exec('generate discarded gene table', self.discarded_gene_table,
                     flag=self.cfg.run.skip_gene_discarded_table, is_sub=True):
            return True

        if self.exec('generate mapping table', self.mapping_table, flag=self.cfg.run.skip_mapping_table, is_sub=True):
            return True

        if self.exec('generate discarded mapping table', self.discarded_mapping_table,
                     flag=self.cfg.run.skip_mapping_discarded_table, is_sub=True):
            return True

        if self.exec('generate ambiguous table', self.ambiguous_table, flag=self.cfg.run.skip_ambiguous_table,
                     is_sub=True):
            return True

        if self.exec('generate raw table', self.raw_table, flag=self.cfg.run.skip_raw_table, is_sub=True):
            return True

        if self.exec('generate regulation table', self.regulation_table, flag=self.cfg.run.skip_regulation_table,
                     is_sub=True):
            return True

        if self.exec('generate annotation tables', self.annotation_tables, flag=self.cfg.run.skip_annotation_tables,
                     is_sub=True):
            return True

        if self.exec('generate mapping .tsv', self.mapping_to_tsv, flag=self.cfg.run.skip_mapping_tsv, is_sub=True):
            return True

        return False

    def fill_table(self, table: str, fields: str, values: list[tuple], indices=None):
        """
        Helper function to create, define and fill a table in the SQLite database. If the table already exists, it and
        its content are first dropped before being created and filled from scratch.
        :param table: table name
        :param fields: comma separated field definitions, e.g. 'identifier TEXT PRIMARY KEY, database_id TEXT'
        :param values: list of rows to be entered into the table, the columns of each row must match the number and
                       order of columns in the fields definition
        :param indices: dictionary of (index_name, is_unique): [column 1, column 2,...]
        """
        # open the database
        connection = sqlite3.connect(self.cfg.paths.database)
        cursor = connection.cursor()
        # drop the table and re-create it from scratch
        cursor.execute(f'DROP TABLE IF EXISTS {table}')
        cursor.execute(f'CREATE TABLE {table} ({fields})')
        # create the specified indices
        if indices:
            for (index, is_unique), columns in indices.items():
                unique = 'UNIQUE ' if is_unique else ''
                cursor.execute(f'CREATE {unique}INDEX {index} ON {table} (' + ', '.join(columns) + ')')
        # enter the provided data into the table
        if values:
            # ?, ?, ?, ?,...
            value_placeholder = ', '.join(['?'] * len(values[0]))
            cursor.executemany(f'INSERT INTO {table} VALUES ({value_placeholder})', values)
        # save the date and close the database
        connection.commit()
        connection.close()

    def gene_table(self):
        """ Builds and fills the table with gene information. """
        table = 'genes'
        fields = ['database_id TEXT PRIMARY KEY', 'type TEXT', 'in_regulation INTEGER',
                  'is_regulator INTEGER', 'is_target TEXT', 'aliases TEXT',
                  'identifier_source TEXT', 'diseases TEXT', 'tissues TEXT', 'processes TEXT', 'DSW INTEGER',
                  'Entrez TEXT']
        values = [(name,
                   gene.type,
                   gene.in_regulation,
                   gene.is_regulator,
                   gene.is_target,
                   ', '.join(sorted(gene.aliases.keys())),
                   gene.sources_str(),
                   '; '.join(sorted(gene.diseases)),
                   '; '.join(sorted(gene.tissues)),
                   '; '.join(sorted(gene.processes)),
                   gene.dsw,
                   ', '.join(sorted(gene.entrez)))
                  for name, gene in self.genes.items()]
        # TODO: generate only the actually required indices
        indices = {(f'genes_type_index', False): ['type'],
                   (f'genes_regulation_index', False): ['in_regulation'],
                   (f'genes_type_regulation_index', False): ['type', 'in_regulation']}
        self.fill_table(table, ', '.join(fields), values, indices)

    def discarded_gene_table(self):
        """ Builds and fills the table with information about discarded genes. """
        table = 'genes_discarded'
        fields = ['database_id TEXT PRIMARY KEY', 'type TEXT', 'aliases TEXT', 'identifier_source TEXT', 'Entrez TEXT']
        values = [(name, gene.type, ', '.join(sorted(gene.aliases.keys())), gene.sources_str(),
                   ', '.join(sorted(gene.entrez))) for name, gene in self.genes_discarded.items()]
        self.fill_table(table, ', '.join(fields), values)

    def mapping_table(self):
        """ Builds and fills the table with mapping information. """
        table = 'mapping'
        fields = ['identifier TEXT PRIMARY KEY', 'database_id TEXT', 'sources TEXT',
                  'FOREIGN KEY (database_id) REFERENCES genes(database_id)']
        values = [(identifier, id_info.name, ', '.join(sorted(id_info.sources)))
                  for identifier, id_info in self.mapping.items()]
        # TODO: generate only the actually required indices
        indices = {(f'mapping_source_index', False): ['sources']}
        self.fill_table(table, ', '.join(fields), values, indices)

    def discarded_mapping_table(self):
        """ Builds and fills the table with discarded mapping information. """
        table = 'mapping_discarded'
        fields = ['identifier TEXT PRIMARY KEY', 'database_id TEXT', 'sources TEXT',
                  'FOREIGN KEY (database_id) REFERENCES genes(database_id)']
        values = [(identifier, id_info.name, ', '.join(sorted(id_info.sources)))
                  for identifier, id_info in self.mapping_discarded.items()]
        self.fill_table(table, ', '.join(fields), values)

    def ambiguous_table(self):
        """
        Builds and fills the table with information of ambiguous identifier mapping that is excluded from the main
        mapping table.
        """
        table = 'mapping_ambiguous'
        fields = ['identifier TEXT PRIMARY KEY', 'database_ids TEXT']
        fmt = '{0} => {1}/{2}/{3}/{4}/{5}/{6} - [{7}]'
        values = [(identifier,
                   ';\n'.join(sorted([fmt.format(id_name, id_info.in_manual_fix, id_info.in_main_cleanup,
                                                 id_info.in_main, id_info.in_miRBase, id_info.in_Entrez,
                                                 id_info.in_regulation, ', '.join(sorted(id_info.sources)))
                                      for id_name, id_info in ids.items()])))
                  for identifier, ids in self.ambiguous.items()]
        self.fill_table(table, ', '.join(fields), values)

    def raw_table(self):
        """
        Builds and fills the table with the raw and unprocessed identifier mapping. This can be useful for development
        but is unnecessary for the final database.
        """
        table = 'mapping_raw'
        fields = ['identifier TEXT PRIMARY KEY', 'database_ids TEXT']
        values = [(identifier,
                   ';\n'.join(sorted(['{0} => {1}/{2}/{3}/{4} - [{5}]'.format(id_name,
                                                                              id_info.in_main,
                                                                              id_info.in_miRBase,
                                                                              id_info.in_Entrez,
                                                                              id_info.in_regulation,
                                                                              ', '.join(sorted(id_info.sources)))
                                      for id_name, id_info in ids.items()])))
                  for identifier, ids in self.mapping_raw.items()]
        self.fill_table(table, ', '.join(fields), values)

    def regulation_table(self):
        """ Builds and fills the table with regulatory information. """
        table = 'regulation'
        fields = ['regulator TEXT', 'target TEXT', 'type TEXT', 'predicted INTEGER', 'experimental INTEGER',
                  'PPI_only INTEGER', 'score INTEGER', 'source TEXT',
                  'FOREIGN KEY (regulator) REFERENCES genes (database_id)',
                  'FOREIGN KEY (target) REFERENCES genes (database_id)',
                  'PRIMARY KEY (regulator, target)']
        values = [(regulator,
                   target,
                   regulation.type,
                   regulation.predicted,
                   regulation.experimental,
                   regulation.ppi_only,
                   regulation.score,
                   ', '.join(sorted(regulation.sources)))
                  for (regulator, target), regulation in self.regulation.items()]
        # TODO: generate only the actually required indices
        indices = {(f'regulation_type_index', False): ['type'],
                   (f'regulation_evidence_index', False): ['predicted', 'experimental'],
                   (f'regulation_type_evidence_index', False): ['type', 'predicted', 'experimental'],
                   (f'regulation_type_evidence_regulator_index', False): ['type', 'predicted', 'experimental',
                                                                          'regulator'],
                   (f'regulation_type_evidence_target_index', False): ['type', 'predicted',  'experimental', 'target'],
                   (f'regulation_general_index', True): ['regulator', 'target', 'type'],
                   (f'regulation_source_index', False): ['source'],
                   (f'regulation_type_evidence_source_index', False): ['type', 'predicted', 'experimental', 'source']}
        self.fill_table(table, ', '.join(fields), values, indices)

    def annotation_tables(self):
        """ Builds and fills a table for each annotation type. """
        for annotation_type, annotation_dict in self.annotations.items():
            table = f'annotations_{annotation_type}'
            fields = ['name TEXT PRIMARY KEY', 'database_ids TEXT']
            values = [(annotation, ', '.join(sorted(genes)))
                      for annotation, genes in sorted(annotation_dict.items(), key=lambda x: x[0].lower())]
            indices = {}
            self.fill_table(table, ', '.join(fields), values, indices)

    def mapping_to_tsv(self):
        """ Generates a tab-separated file with identifier mapping. """
        with open(self.cfg.paths.mapping_tsv, 'w') as file:
            file.write('input.id\tdatabase.id\n')
            for identifier, symbol in sorted(self.mapping.items()):
                file.write(f'{identifier}\t{symbol.name}\n')

    # INFORMATION & STATISTICS OUTPUT ----------------------------------------------------------------------------------
    def report(self):
        """
        Outputs the log dictionary in .yaml format. The log dictionary is organised as follows:
        {category: {sub-category: {key: value}}}
        The values can be lists, sets, dictionaries and numbers. By default, only the length of lists, sets and
        dictionaries is reported.
        """
        # List of tuples [(category, sub-category, key),...] of log entries where the elements of the set or list
        # should be reported in addition to its length. For dictionaries the elements are the keys.
        verbose = self.cfg.run.log_verbose                                  # type: set[tuple[str, str, str]]

        with open(self.cfg.paths.log, 'w', encoding='utf-8') as file:
            # iterate over the categories in the order they were entered into the log
            for category, subcategories in self.log_dict.items():
                file.write(f'{category}:\n')
                # iterate over the sub-categories in the order they were entered into the category
                for sub, sub_dict in subcategories.items():
                    file.write(f'\t{sub}:\n')
                    # build the line format to align the values
                    key_len = max(len(k) for k in sub_dict.keys()) + 1
                    int_len = max(len(str(x)) for x in [v if isinstance(v, int) else len(v) for v in sub_dict.values()])
                    if int_len > 6:
                        int_len += 3
                    elif int_len > 3:
                        int_len += 2
                    f = '\t\t{0:<' + str(key_len) + '} {1:>' + str(int_len) + ',}'
                    # iterate over the key, value pairs in alphabetical order of the keys
                    for k, v in sorted(sub_dict.items()):
                        # case: the value is a set or list => by default output the number of elements,
                        # otherwise sorted elements
                        if isinstance(v, set) or isinstance(v, list):
                            if verbose and (category, sub, k) in verbose:
                                file.write((f + ' => {2}\n').format(k + ':', len(v), ', '.join(sorted(v))))
                            else:
                                file.write((f + '\n').format(k + ':', len(v)))
                        # case: the value is a dictionary => by default output the number of keys, otherwise sorted keys
                        elif isinstance(v, dict):
                            if verbose and (category, sub, k) in verbose:
                                file.write((f + ' => {2}\n').format(k + ':', len(v), ', '.join(sorted(v.keys()))))
                            else:
                                file.write((f + '\n').format(k + ':', len(v)))
                        # case: the value is already number
                        else:
                            file.write((f + '\n').format(k + ':', v))

    def get_stats(self) -> dict[str, dict[str, int | list[str]]]:
        """
        :return: dictionary with some statistical information about the database
        """
        # initialise the stats dictionary
        stats = {'interactions': defaultdict(int), 'identifiers': defaultdict(int),
                 'annotations': defaultdict(int), 'update date': str(date.today())}
        # total number of interactions in the database
        stats['interactions']['total'] = len(self.regulation)
        # compile the set of all regulation data sources used in the database
        sources = set()
        for reg in self.regulation.values():
            stats['interactions'][reg.type] += 1
            sources.update(reg.sources)
        stats['interactions']['sources'] = sorted(sources)
        # total number of identifiers and symbols contained in the mapping, not counting discarded mappings
        stats['identifiers']['recognized identifiers'] = len(self.mapping)
        stats['identifiers']['symbols, total'] = len(self.genes)
        # compile the ID systems for the identifier mapping
        systems = set()
        for symbol in self.genes.values():
            stats['identifiers']['symbols, ' + symbol.type] += 1
            systems.update(symbol.sources_str(['main merge', 'main fix', 'manual fix', 'in-house processing']).split(', '))
        stats['identifiers']['ID systems'] = sorted(systems)
        # get the sources underpping the identifier mapping
        sources = set()
        if not self.cfg.run.skip:
            sources.add(self.cfg.paths.main.label())
        if not self.cfg.run.skip_miRBase:
            sources.add(self.cfg.paths.miRBase.label())
        if not self.cfg.run.skip_Entrez:
            sources.add(self.cfg.paths.Entrez.label())
        stats['identifiers']['sources'] = sorted(sources)
        # count the annotations by type
        for annotation_type, annotations in self.annotations.items():
            stats['annotations'][annotation_type] = len(annotations)
        # transform the default dictionaries to normal dictionaries for later outputting with YAML
        for key in ['interactions', 'annotations', 'identifiers']:
            stats[key] = dict(stats[key])

        return stats
