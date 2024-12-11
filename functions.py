import os
import re
import datetime
from time import time
from typing import Dict, Iterable, List, Set, Tuple


class LongestOverlap:
    def __init__(self, identifier: str, symbols: Iterable[str], suffix=False, try_case_insensitive=True):
        """
        Computes the longest overlap between an identifier and each symbol in the input list.
        :param identifier: a gene identifier
        :param symbols: lists of gene symbols
        :param suffix: True if the overlap should be the longest common suffix, False for the longest common prefix
        :param try_case_insensitive: True if the overlap should ignore the case, False otherwise
        """
        self.type = 'suffix' if suffix else 'prefix'        # type: str
        self.identifier = identifier                        # type: str
        # longest common prefix (or suffix) and its length
        self.length = 0                                     # type: int
        self.overlap = ''                                   # type: str
        # symbols that have the longest common prefix (or suffix) in common with the identifier
        self.matches = []                                   # type: List[str]

        # nothing to do if there are no symbols given
        if not symbols:
            return

        # reverse the identifier and the symbols to compute the suffix (however, do not change the versions saved in
        # the class itself)
        if suffix:
            identifier = identifier[::-1]
            symbols = [x[::-1] for x in symbols]

        # identify the longest common prefix (or suffix) between the identifier and each input symbol. if
        # case-insensitive search is enabled, repeat the procedure with identifier and symbols cast to upper case to
        # see if that improves the result
        for to_process_list in [symbols, [x.upper() for x in symbols]] if try_case_insensitive else [symbols]:
            # a unique match was found in the previous (case-sensitive) iteration => best result
            if len(self.matches) == 1:
                break

            # compute the overlap between the identifier and each symbol
            overlap = {(symbol, os.path.commonprefix([identifier.upper() if to_process_list else identifier,
                                                      symbol.upper() if to_process_list else symbol]))
                       for symbol in symbols}
            # sort the overlaps by their length so that the longest ones are first
            sorted_overlap = sorted({(len(common), common, symbol) for symbol, common in overlap}, reverse=True)

            # only save the common prefix (or suffix) if it is longer than the previously saved one
            if sorted_overlap[0][0] < self.length:
                continue
            self.length = sorted_overlap[0][0]
            self.overlap = sorted_overlap[0][1]
            # save all symbols that share the longest overlap with the identifier
            self.matches = list({symbol for length, _, symbol in sorted_overlap if length == self.length})

            # since the identifier and symbols were reversed to compute the longest common suffix, the overlap and
            # matches need to be reversed back into the original direction before saving them
            if suffix:
                self.overlap = self.overlap[::-1]
                self.matches = [symbol[::-1] for symbol in self.matches]


def to_upper(match):
    """ Converts a regex match to upper case. """
    return match.group(0).upper()


def time_stamp(start):
    """
    :param start: start time (generated with time())
    :return: formatted time difference between the start time and the current time
    """
    h, m, s = str(datetime.timedelta(seconds=time() - start)).split(':')
    stamp = [f'{int(h)}h'] if int(h) else []
    if int(m):
        stamp.append(f'{int(m)}min')
    if float(s):
        stamp.append(f'{round(float(s), 1)}s')

    return ' '.join(stamp) if stamp else '0s'


def is_miRBase_accession(identifier: str) -> bool:
    """
    :return: True if the identifier is a miRBase accession (e.g. MI0000060), False otherwise
    """
    return re.fullmatch('(MI)[0-9]{7}', identifier.upper().strip()) is not None


def is_miRBase_ID(identifier: str, species_tag: str) -> bool:
    """
    Checks if the identifier is a miRBase identifier, e.g. hsa-mir-123.
    :param identifier: the identifier to check
    :param species_tag: e.g. 'hsa' for humans, 'mmu' for mouse
    :return: True if the identifier is a miRBase identifier, False otherwise
    """
    species_tag = species_tag.lower()
    # identifiers can have inconsistent usage of capital letters and no whitespace
    identifier = identifier.strip().lower().replace(' ', '')
    # the identifiers have either 'mir' or 'let' in them
    if 'mir' not in identifier and 'let' not in identifier:
        return False
    # some databases have a typo in 'hsa' (species tag for humans)
    if species_tag == 'hsa':
        identifier = identifier.replace('has', 'hsa')

    mir_type = 'mir' if 'mir' in identifier else 'let'

    # case: hsa-mir123 instead of hsa-mir-123
    if mir_type + '-' not in identifier:
        identifier = identifier.replace(mir_type, mir_type + '-')
    # TODO: case: wrong species, fix later
    elif not identifier.startswith(species_tag) and not identifier.startswith(mir_type):
        components = identifier.split('-')
        identifier = species_tag + '-' + '-'.join(components[1:])
    # remove trailing -
    identifier = identifier.strip('-')

    return re.fullmatch(f'({species_tag}-)?({mir_type})(-[a-z0-9*]+)(-[0-9a-z*]+)?(-[0-9a-z*]+)?(-[0-9]+)?',
                        identifier) is not None


def miRBase_ID_to_symbol(identifier: str, species_tag: str) -> Dict[str, str]:
    """
    Converts a miRBase identifier (e.g. hsa-mir-123) to a miRNA symbol (e.g. MIR123)
    :param identifier: identifier to convert
    :param species_tag: e.g. 'hsa' for humans, 'mmu' for mouse
    :return: dictionary of conversion candidates (key) with corresponding conversion method (values)
    """
    # pre-process the identifier, remove the species tag and cast it to upper case (default for miRNA symbols)
    identifier, _ = process_mirBase_ID(identifier, species_tag)
    # not actually a miRBase ID
    if not identifier:
        return dict()
    identifier = identifier.replace(species_tag + '-', '').upper()
    # split identifier into its components and pre-process them into prefix and suffix
    parts = identifier.split('-')
    prefix = 'MIRLET' if parts[0] == 'LET' else 'MIR'
    suffix = '-'.join(parts[1:])
    # go over different ways of converting the miRBase ID depending on the structure of the identifier and return
    # the correspondingly generated miRNA symbols
    if re.fullmatch('[0-9]+', suffix) is not None:
        return {prefix + suffix: '1234',
                prefix + suffix + '-1': '1234, added -1',
                prefix + suffix + 'A1': '1234, added A1',
                prefix + suffix + 'A': '1234, added A'}
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if re.fullmatch('[0-9]+-[0-9]+', suffix) is not None:
        return {prefix + suffix: '1234-1',
                prefix + parts[1] + alphabet[int(parts[2]) - 1]: '1234-1, added A',
                prefix + parts[1] + alphabet[int(parts[2]) - 1] + '-1': '1234-1, added A-1',
                prefix + parts[1] + alphabet[int(parts[2]) - 1] + 'A1': '1234-1, added A1',
                prefix + parts[1]: '1234-1, removed -1'}
    if re.fullmatch('[0-9]+-[0-9]+-[0-9]+', suffix) is not None:
        return {prefix + suffix: '1234-1-201'}
    if re.fullmatch('[0-9]+[A-Z]+', suffix) is not None:
        return {prefix + suffix: '1234A',
                prefix + suffix + '-1': '1234A, added -1',
                prefix + suffix + '1': '1234A, added 1'}
    if re.fullmatch('[0-9]+[A-Z]+-[0-9]+', suffix) is not None:
        return {prefix + suffix: '1234A-1',
                prefix + suffix.replace('-', ''): '1234A-1, removed -'}
    if re.fullmatch('[0-9]+[A-Z]+[0-9]+', suffix) is not None:
        return {prefix + suffix: '1234A1'}
    return {prefix + suffix: 'remaining'}


def process_mirBase_ID(identifier: str, species_tag: str) -> Tuple[str, Set[str]]:
    """
    Processes a miRBase ID (e.g. hsa-mir-123) to account for common spelling mistakes, variances in spelling or
    upper or lower case spelling.
    :param identifier: identifier to process
    :param species_tag: e.g. 'hsa' for humans, 'mmu' for mouse
    :return: the processed identifier and a set of spelling versions of the identifier
    """
    # only process miRBase IDs
    if not is_miRBase_ID(identifier, species_tag):
        return '', set()
    # initialise with the input version of the identifier (without potential surrounding white space)
    identifier = identifier.strip()
    transformed = {identifier, identifier.lower()}
    # remove whitespace
    identifier = identifier.replace(' ', '')
    transformed.add(identifier)
    # case to lower case for further case independent processing and because miRBase IDs should be lowercase anyway
    identifier = identifier.lower()
    # fix a common spelling mistake in the human species tag
    if species_tag == 'hsa' and identifier.startswith('has'):
        identifier = identifier.replace('has', 'hsa')
        transformed.add(identifier)
    # wrong species
    mir_type = 'mir' if 'mir' in identifier else 'let'
    if not identifier.startswith(species_tag) and not identifier.startswith(mir_type):
        return '', set()
    # add the species tag if it's missing
    if not identifier.startswith(species_tag):
        identifier = species_tag + '-' + identifier
        transformed.add(identifier)
    # case: mir1231 or hsa-mir223
    if mir_type + '-' not in identifier:
        identifier = identifier.replace(mir_type, mir_type + '-')
        transformed.add(identifier)
    # create a version without the species tag
    identifier_without = identifier.replace(species_tag + '-', '')
    transformed.add(identifier_without)
    # for both versions go over various combinations of spellings, components and additional modifiers and add them
    # to the set of possible versions for the identifier
    for i in [identifier, identifier_without]:
        for old_tag, new_tag in [('', ''), (species_tag, species_tag.capitalize()), (species_tag, species_tag.upper())]:
            x = i.replace(old_tag, new_tag)
            for old, new in [('', ''), ('mir', 'miR'), ('mir', 'MiR'), ('mir', 'MIR'), ('let', 'LET')]:
                y = x.replace(old, new)
                transformed.add(y)
                for strip in ['*', '-3p', '-5p']:
                    z = y[:-len(strip)] if y.endswith(strip) else y
                    transformed.add(z)
    # remove tags from the identifier
    for strip in ['*', '-3p', '-5p']:
        identifier = identifier[:-len(strip)] if identifier.endswith(strip) else identifier
        transformed.add(identifier)

    return identifier.lower(), transformed


def is_mirna_symbol(identifier: str) -> bool:
    """ Checks if the identifier is a miRBase symbol """
    return re.fullmatch('(MIR)(N)?(LET)?[0-9]+(HG)?-?[0-9A-Z]*-?[0-9A-Z*]*-?[0-9]*', identifier.upper()) is not None


def is_invalid(identifier: str) -> str:
    """
    :return: error string if there is some basic issue with the identifier, empty string if the identifier seems fine
    """
    identifier = identifier.lower().strip()
    if identifier in ['n/a', 'na', 'unknown', 'none', 'n_a', 'n.a.', '']:
        return 'empty identifier'
    for x in ['https', 'unknown', 'http:']:
        if x in identifier:
            return 'invalid component'
    if identifier == 'eco:' or (identifier.startswith('eco:') and identifier[4:].split('|')[0].isdecimal()):
        return 'ECO evidence'
    return ''
