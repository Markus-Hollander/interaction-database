import re
import os
import yaml
from datetime import date
from typing import Dict, List, Set, Tuple, Union
            

class Formatting:
    """ Holds information for fixing and unifying the spelling and formatting of gene annotations. """
    def __init__(self, cfg: Dict):
        self.invalid = set(cfg.get('invalid', []))                                  # type: Set[str]
        # removal or replacement
        self.remove = cfg.get('remove', [])                                         # type: List[str]
        self.replace_single = dict()                                                # type: Dict[str, re.Pattern]
        self.replace_multiple = cfg.get('replace, multiple', dict())                # type: Dict[str, str]
        self.replace_regex = dict()                                                 # type: Dict[str, re.Pattern]
        # case adjustment
        self.upper_complex = []                                                     # type: List[re.Pattern]
        self.upper_simple = dict()                                                  # type: Dict[str, re.Pattern]
        self.lower = dict()                                                         # type: Dict[str, re.Pattern]
        self.capitalize = dict()                                                    # type: Dict[str, re.Pattern]
        self.mixed = dict()                                                         # type: Dict[str, re.Pattern]

        # pre-compile all REGEX patterns so that it only has to be done once
        for pat, repl in cfg.get('replace, single', dict()).items():
            regex = '(?<![0-9a-zA-Z])(' + ' '.join(pat.strip().split()).lower().replace('(', '\(').replace(')', '\)')
            regex += ')(?![a-zA-Z])'
            pattern = re.compile(regex, flags=re.IGNORECASE)
            self.replace_single[repl] = pattern

        for pattern, repl in cfg.get('replace, regex', dict()).items():
            self.replace_single[repl] = re.compile(pattern, flags=re.IGNORECASE)

        for low in set(cfg.get('lower, single', [])):
            regex = '(?<![0-9a-zA-Z])(' + ' '.join(low.strip().split()).lower() + ')(?![a-zA-Z])'
            pattern = re.compile(regex, flags=re.IGNORECASE)
            self.lower[low.lower()] = pattern

        for up in set(cfg.get('upper, complex', [])):
            self.upper_complex.append(re.compile(up, flags=re.IGNORECASE))

        for up in set(cfg.get('upper, simple', [])):
            regex = '(?<![0-9a-zA-Z])(' + ' '.join(up.strip().split()).lower() + ')(?![a-zA-Z])'
            pattern = re.compile(regex, flags=re.IGNORECASE)
            self.upper_simple[up.upper()] = pattern

        for cap in set(cfg.get('capitalize, single', [])):
            regex = '(?<![0-9a-zA-Z])(' + ' '.join(cap.strip().split()).lower() + ')(?![a-zA-Z])'
            pattern = re.compile(regex, flags=re.IGNORECASE)
            self.capitalize[cap.capitalize()] = pattern

        for sep in ['+', '-']:
            for cap in set(cfg.get(f'capitalize, {sep}', [])):
                cap = cap.strip().lower()
                regex = '(?<![0-9a-zA-Z])(' + f'[\\s{sep}]'.join(cap.split()).lower() + ')(?![a-zA-Z])'
                pattern = re.compile(regex, flags=re.IGNORECASE)
                self.capitalize[sep.join([x.capitalize() for x in cap.split()])] = pattern

        for mixed in set(cfg.get('mixed', [])):
            regex = '(?<![0-9a-zA-Z])(' + ' '.join(mixed.strip().split()).lower() + ')(?![a-zA-Z])'
            pattern = re.compile(regex, flags=re.IGNORECASE)
            self.mixed[mixed] = pattern


class File:
    """ Information about the layout of a single general raw data input file and how to automatically process it. """
    def __init__(self, file_path: str, cfg: Dict[str, str]):
        self.name = cfg.get('name', '')                                             # type: str
        self.file_name = cfg.get('file name', '')                                   # type: str
        self.path = file_path                                                       # type: str
        self.version = str(cfg.get('version', ''))                                  # type: str
        # LAYOUT
        self.delimiter = cfg.get('delimiter', '\t')                                 # type: str
        self.header = cfg.get('header', True)                                       # type: bool
        self.remaining_cols = cfg.get('remaining cols', False)                      # type: bool
        self.min_cols = cfg.get('minimum columns', 0)                               # type: int

    def label(self) -> str:
        """ Dynamically generate the label. """
        if not self.name:
            return '.'.join(self.file_name.split('.')[:-1]).replace('_', ' ')
        label = self.name
        if self.version:
            label += f' ({self.version})'
        return label


class Annotation(File):
    """ Information about the layout of a single raw annotation data input file and how to automatically process it. """
    def __init__(self, file_path: str, cfg: Dict):
        super().__init__(file_path, cfg)
        self.type = cfg.get('type', 'MISSING')                                      # type: str
        self.annotation = cfg.get('annotation', 'MISSING')                          # type: str
        self.gene = cfg.get('gene', 'MISSING')                                      # type: str
        self.split = cfg.get('annotation split', '')                                # type: str


class Regulation(File):
    """ Information about the layout of a single raw regulation data input file and how to automatically process it. """
    def __init__(self, source: str, file_path: str, cfg: Dict):
        super().__init__(file_path, cfg)
        self.source = source                                                        # type: str
        self.evidence = cfg.get('evidence', 'MISSING')                              # type: str
        self.regulator_cols = cfg.get('regulator columns', [])                      # type: List[str]
        self.regulator_type = cfg.get('regulator type', '')                         # type: str
        self.target_cols = cfg.get('target columns', [])                            # type: List[str]
        self.target_type = cfg.get('target type', '')                               # type: str
        self.score_col = cfg.get('score column', '')                                # type: str
        self.score_divisor = cfg.get('score divisor', 1.0)                          # type: float
        self.species_cols = cfg.get('species columns', [])                          # type: List[str]
        self.ppi_only = cfg.get('ppi only', False)                                  # type: bool


class PathConfig:
    """ File and directory paths for retrieving and storing data. """
    def __init__(self):
        # MAPPING
        self.mapping_dir = ''               # type: str
        self.mapping_merge_dir = ''         # type: str
        self.main = None                    # type: Union[File, None]
        self.miRBase = None                 # type: Union[File, None]
        self.Entrez = None                  # type: Union[File, None]
        self.fixes = None                   # type: Union[File, None]
        # ADDITIONAL DATA
        self.gene_type = None               # type: Union[File, None]
        self.DSW = None                     # type: Union[File, None]
        self.annotations = []               # type: List[Annotation]
        self.regulation = []                # type: List[Regulation]
        # OUTPUT
        self.log = ''                       # type: str
        self.database = ''                  # type: str
        self.mapping_tsv = ''               # type: str


class RunConfig:
    """ Configuration for enabling or disabling certain parts of the database generation script. """
    def __init__(self):
        self.skip = ''                          # type: str
        # MAPPING
        self.skip_merge = ''                    # type: str
        self.skip_fixes = ''                    # type: str
        self.skip_miRBase = ''                  # type: str
        self.skip_Entrez = ''                   # type: str
        self.skip_untangle = ''                 # type: str
        self.skip_final_cleanup = ''            # type: str
        # MAPPING - CONFIG
        self.main_cols = set()                  # type: Set[str]
        self.Entrez_cols = set()                # type: Set[str]
        # ADDITIONAL DATA
        self.skip_gene_type = ''                # type: str
        self.skip_DSW = ''                      # type: str
        self.skip_annotations = ''              # type: str
        # REGULATION
        self.skip_regulation = ''               # type: str
        self.skip_regulation_untangle = ''      # type: str
        self.skip_regulation_resolve = ''       # type: str
        # DATABASE
        self.skip_database = ''                 # type: str
        self.skip_gene_table = ''               # type: str
        self.skip_gene_discarded_table = ''     # type: str
        self.skip_mapping_table = ''            # type: str
        self.skip_mapping_discarded_table = ''  # type: str
        self.skip_ambiguous_table = ''          # type: str
        self.skip_raw_table = ''                # type: str
        self.skip_regulation_table = ''         # type: str
        self.skip_annotation_tables = ''        # type: str
        self.skip_mapping_tsv = ''              # type: str
        # REPORT
        self.skip_log = ''                      # type: str
        self.log_verbose = set()                # type: Set[Tuple[str, str, str]]


class JobConfig:
    """ Configuration and information for generating a single database version of a species, e.g. human version 1. """
    def __init__(self, version: str, species: str, species_cfg: Dict[str, str], format_cfg: Formatting):
        self.version = version                                  # type: str
        # GENERAL SPECIES INFO
        self.species = species                                  # type: str
        self.species_full = species_cfg['name']                 # type: str
        self.miRNA_prefix = species_cfg['miRNA prefix']         # type: str
        self.taxon_id = species_cfg['taxon ID']                 # type: str
        # TERMS
        self.no_symbol = 'NO_SYMBOL'                            # type: str
        # FORMATTING
        self.format = format_cfg                                # type: Formatting
        # PATHS and RUN FLAGS
        self.run = RunConfig()                                  # type: RunConfig
        self.paths = PathConfig()                               # type: PathConfig


class Config:
    """ Stores and manages all configurations for central and easy access with auto-complete. """
    def __init__(self, run_config_path: str, db_config_path: str, annotation_config_path: str):
        # CONFIGS
        with open(run_config_path, 'r') as file:
            run_config = yaml.safe_load(file)                                                   # type: Dict
        with open(db_config_path, 'r') as file:
            db_config = yaml.safe_load(file)                                                    # type: Dict
        with open(annotation_config_path, 'r') as file:
            annotation_config = yaml.safe_load(file)                                            # type: Dict

        # DIRECTORIES
        self.directories = db_config['directories']                                             # type: Dict[str, str]
        self.databases_dir = db_config['directories']['database']                               # type: str
        self.processing_dir = db_config['directories']['processing']                            # type: str
        # create the database directory if it does not exist yet, the others are created later
        os.makedirs(self.databases_dir, exist_ok=True)
        # path to information file about all created databases
        self.report_file = os.path.join(self.databases_dir, 'database_overview.yaml')           # type: str

        # FLAGS
        self.delete_database = run_config['delete old database']                                # type: bool

        # FORMATTING
        format_cfg = Formatting(annotation_config)                                              # type: Formatting

        # JOBS
        self.jobs = []                                                                          # type: List[JobConfig]

        # generate the configurations for each job (defined by version and species)
        # v: version
        for v, version_dict in sorted(db_config['versions'].items()):
            # s: species
            for s, config in sorted(version_dict.items()):
                # initialise the job specific config with species information
                job = JobConfig(v, s, db_config['species'][s], format_cfg)
                # OUTPUT PATHS
                # log path
                log_dir = self.job_path(v, s, 'logs')
                today = str(date.today())
                existing_logs = sum(1 for directory in os.listdir(log_dir) if directory.startswith(today))
                job.paths.log = os.path.join(log_dir, f'{today}_log_{existing_logs}.yaml')
                # log configuration
                job.run.log_verbose = run_config['report']['verbose']
                # database
                job.paths.database = os.path.join(self.databases_dir, f'{s}_version_{v}.db')
                # additional output files
                job.paths.mapping_tsv = os.path.join(self.job_path(v, s, 'mapping'),
                                                     f'generated_mapping_{s}_version{v}.tsv')

                # MAPPING
                if 'mapping' in config:
                    job.paths.mapping_dir = self.job_path(v, s, 'mapping')
                    job.paths.mapping_merge_dir = self.job_path(v, s, 'mapping merge')

                    if not os.listdir(job.paths.mapping_merge_dir):
                        job.run.skip_merge = 'no files to merge'

                    cfg = config['mapping']
                    job.paths.main, job.run.skip = self.job_path(v, s, 'mapping', cfg, 'main')
                    job.paths.fixes, job.run.skip_fixes = self.job_path(v, s, 'mapping', cfg, 'fixes')
                    job.paths.miRBase, job.run.skip_miRBase = self.job_path(v, s, 'mapping', cfg, 'miRBase')
                    job.paths.Entrez, job.run.skip_Entrez = self.job_path(v, s, 'mapping', cfg, 'Entrez')

                    if 'columns' not in cfg['main'] or not cfg['main']['columns']:
                        job.run.skip = 'no mapping columns'
                    else:
                        job.run.main_cols = set(cfg['main']['columns'])

                    if 'Entrez columns' in cfg['main'] and cfg['main']['Entrez columns']:
                        job.run.Entrez_cols = set(cfg['main']['Entrez columns'])
                else:
                    job.run.skip = 'no mapping configuration'

                # ANNOTATIONS
                if 'additional' in config:
                    cfg = config['additional']
                    job.paths.gene_type, job.run.skip_gene_type = self.job_path(v, s, 'annotations', cfg, 'gene types')
                    job.paths.DSW, job.run.skip_DSW = self.job_path(v, s, 'annotations', cfg, 'DSW')
                else:
                    job.run.skip_gene_type = 'no gene type configuration'
                    job.run.skip_DSW = 'no DSW configuration'

                if 'annotations' in config:
                    for file_name, cfg in config['annotations'].items():
                        file_path, error = self.job_path(v, s, 'annotations', key=file_name)
                        anno = Annotation(file_path, cfg)
                        if error or 'MISSING' in [anno.annotation, anno.gene]:
                            continue
                        if anno.type not in ['tissue', 'disease', 'process']:
                            continue
                        job.paths.annotations.append(anno)
                if not job.paths.annotations:
                    job.run.skip_annotations = 'no annotations'

                # REGULATION
                if 'regulation' in config:
                    for source, cfg in config['regulation'].items():

                        reg_file, error = self.job_path(v, s, 'regulation', config['regulation'], source)
                        reg = Regulation(source, reg_file.path, cfg)
                        if error or reg.evidence == 'MISSING' or not reg.regulator_cols or not reg.target_cols:
                            continue
                        job.paths.regulation.append(reg)
                if not job.paths.regulation:
                    job.run.skip_regulation = 'no regulatory information'

                # RUN FLAGS
                if v in run_config['skip versions']:
                    job.run.skip = 'skip version'
                elif s in run_config['skip species']:
                    job.run.skip = 'skip species'
                elif v in run_config['skip specific'] and s in run_config['skip specific'][v]:
                    job.run.skip = 'skip specific'
                # mapping flags
                if not run_config['mapping']['merge main']:
                    job.run.skip_merge = 'disabled'
                if not run_config['mapping']['load fixes']:
                    job.run.skip_fixes = 'disabled'
                if not run_config['mapping']['load additional Entrez']:
                    job.run.skip_Entrez = 'disabled'
                if not run_config['mapping']['load additional miRBase']:
                    job.run.skip_miRBase = 'disabled'
                if not run_config['mapping']['untangle']:
                    job.run.skip_untangle = 'disabled'
                if not run_config['mapping']['final cleanup']:
                    job.run.skip_final_cleanup = 'disabled'
                # annotation flags
                if not run_config['annotations']['load gene types']:
                    job.run.skip_gene_type = 'disabled'
                if not run_config['annotations']['load DSW']:
                    job.run.skip_DSW = 'disabled'
                if not run_config['annotations']['load annotations']:
                    job.run.skip_annotations = 'disabled'
                # regulation flags
                if not run_config['regulation']['load']:
                    job.run.skip_regulation = 'disabled'
                    job.run.skip_regulation_resolve = 'disabled'
                    job.run.skip_regulation_untangle = 'disabled'
                if not run_config['regulation']['untangle']:
                    job.run.skip_regulation_untangle = 'disabled'
                if not run_config['regulation']['resolve']:
                    job.run.skip_regulation_resolve = 'disabled'
                # output flags
                if not run_config['database']['output']:
                    job.run.skip_database = 'disabled'
                if not run_config['database']['gene table']:
                    job.run.skip_gene_table = 'disabled'
                if not run_config['database']['discarded gene table']:
                    job.run.skip_gene_discarded_table = 'disabled'
                if not run_config['database']['mapping table']:
                    job.run.skip_mapping_table = 'disabled'
                if not run_config['database']['discarded mapping table']:
                    job.run.skip_mapping_discarded_table = 'disabled'
                if not run_config['database']['ambiguous table']:
                    job.run.skip_ambiguous_table = 'disabled'
                if not run_config['database']['raw table']:
                    job.run.skip_raw_table = 'disabled'
                if not run_config['database']['regulation table']:
                    job.run.skip_regulation_table = 'disabled'
                if not run_config['database']['annotation tables']:
                    job.run.skip_annotation_tables = 'disabled'
                if not run_config['database']['mapping tsv']:
                    job.run.skip_mapping_tsv = 'disabled'
                if not run_config['report']['output']:
                    job.run.skip_log = 'disabled'

                self.jobs.append(job)

        # get lists of all configured versions and database paths
        self.versions = {job.version for job in self.jobs if not job.run.skip}                      # type: Set[str]
        self.database_paths = {job.paths.database for job in self.jobs if not job.run.skip}         # type: Set[str]

    def job_path(self, version: str, species: str, dir_key: str, cfg=None,
                 key='') -> Union[str, Tuple[str, str], Tuple[File, str]]:
        """
        Generate a file or directory path in the processing directory for a job, defined by version and species.
        :param version: database version, becomes part of the directory path
        :param species: species name, becomes part of the directory path
        :param dir_key: key to retrieve the sub-category name from the directory name dictionary
        :param cfg: configuration dictionary
        :param key: key in cfg if cfg is given, or file name if cfg is not given
        :return:
            case 1: directory path if cfg and key are not given
            case 2: (file path, error string) if cfg and/or key are given
        """
        # generate the directory path: processing_dir/version_X/species/sub-directory/
        directory = os.path.join(self.processing_dir, f'version_{version}', species, self.directories[dir_key])
        # create the directory if it does not exist yet
        os.makedirs(directory, exist_ok=True)
        # only return the directory path if no file information is provided
        if not cfg and not key:
            return directory
        # otherwise generate the file path and check if it points to an existing file
        # case: a config dictionary is given but the key is not in the config dictionary => error
        if cfg and key not in cfg:
            return '', 'not in configuration'
        if cfg and 'file name' not in cfg[key]:
            return '', 'file name missing'
        # if only the key is given, the key is actually the desired file name, otherwise the file name is in the
        # config dictionary
        file_path = os.path.join(directory, cfg[key]['file name'] if cfg else key)
        try:
            # case: no error occurred => return the file path and an empty error string
            with open(file_path, 'r') as _:
                return File(file_path, cfg[key]) if cfg else file_path, ''
        # case: the file path does not point to an existing file => error
        except IOError:
            return File(file_path, cfg[key]) if cfg else file_path, f'{key} file does not exist'
