from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional

from functions import is_miRBase_accession


@dataclass(slots=True)
class Gene:
    """ Represents a gene and information associated with the gene. """
    # name of the gene, typically the gene symbol
    name: str = ''
    # the gene type is one of transcription factor, miRNA or gene (default)
    type: str = 'gene'
    # True if there is at least one interaction involving the gene, False otherwise
    in_regulation: bool = False
    is_regulator: bool = False
    is_target: bool = False
    # all identifiers that are mapped to this gene, as well as the sources underpinning the mappings
    # {identifier: {sources of identifier -> gene mapping},...}
    aliases: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    # set of Entrez IDs associated with this gene, used for some TFmiR functionalities
    # (only filled in at the end, for the raw mapping use get_Entrez_ids)
    entrez: set[str] = field(default_factory=set)
    # sets of tissue, biological process and disease annotations associated with the gene
    tissues: set[str] = field(default_factory=set)
    processes: set[str] = field(default_factory=set)
    diseases: set[str] = field(default_factory=set)
    # disease spectrum width (only relevant for miRNAs)
    dsw: float = 0.0

    def get_miRBase_accessions(self) -> set[str]:
        """ :return: set of all miRBase accessions associated with the gene """
        return {a for a in self.aliases.keys() if is_miRBase_accession(a)}

    def get_Entrez_ids(self, ignore=None) -> set[str]:
        """
        :param ignore: set of identifiers that should be excluded from the list of Entrez ids
        :return: set of all Entrez IDs associated with the gene
        """
        if ignore:
            return {a for a in self.aliases.keys() if a.isdigit() and a not in ignore}
        return {a for a in self.aliases.keys() if a.isdigit()}

    def aliases_str(self, details=False) -> str:
        """
        :param details: True if the sources of the alias mapping should be included in the string, False otherwise
        :return: formatted string containing all aliases associated with the gene
        """
        if details:
            return ';\n'.join(sorted(['{0} ({1})'.format(k,
                                                         ', '.join(v)) for k, v in self.aliases.items()]))
        return ', '.join(sorted(self.aliases.keys()))

    def sources_str(self, exclude=None) -> str:
        """
        :param exclude: set of sources that should not be included in the formatted string
        :return: formatted string of all sources mentioning the gene
        """
        sources = set()
        for mapping_sources in self.aliases.values():
            sources.update(mapping_sources)
        for to_remove in exclude if exclude else []:
            sources.discard(to_remove)
        sources = {s.replace(' (formerly Entrezgene) ',
                             '') for s in sources}
        return ', '.join(sorted(sources))


@dataclass(slots=True)
class Symbol:
    """ Represents a gene symbol in the raw or final identifier -> gene symbol mapping. """
    # INFO
    name: str = ''
    type: str = 'gene'
    sources: Optional[set[str]] = None
    # SOURCE FLAGS
    in_manual_fix: bool = False
    in_main_cleanup: bool = False
    in_miRBase: bool = False
    in_main: bool = False
    in_Entrez: bool = False
    in_regulation: bool = False

    def __post_init__(self):
        if not self.sources:
            self.sources = set()


@dataclass(slots=True)
class Regulation:
    """ Stores information for a single interaction. """
    type: str = ''
    experimental: bool = False
    predicted: bool = False
    sources: set[str] = field(default_factory=set)
    score: float = 1.0
    ppi_only: bool = False
