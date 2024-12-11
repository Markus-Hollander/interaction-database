from collections import defaultdict
import csv
import os
import sys
from time import time
import yaml

from config import Config
from database import DataBaseBuilder
from functions import time_stamp


def generate():
    """ Runs all required database building steps for the specified configurations. """
    config = Config('config/run.yaml', 'config/database.yaml', 'config/annotation_formatting.yaml')     # type: Config

    # if specified, delete old databases that are changed in any of the valid jobs
    if config.delete_database:
        for p in config.database_paths:
            if os.path.isfile(p):
                os.remove(p)

    # SETUP
    # decrease the maxInt value by factor 10
    # as long as the OverflowError occurs.
    maxInt = sys.maxsize
    while True:
        try:
            csv.field_size_limit(maxInt)
            break
        except OverflowError:
            maxInt = int(maxInt / 10)

    # STAT SETUP
    stats = defaultdict(lambda: defaultdict(dict))          # type: dict[str, dict[str, dict]]
    # load information from previous runs to not lose it (if the file already exists from a previous run)
    try:
        with open(config.report_file, 'r') as stat_file:
            previous_stats = yaml.safe_load(stat_file)
        if previous_stats:
            for species, species_info in previous_stats.items():
                for version, version_info in species_info.items():
                    stats[species][version] = version_info
    # case: file does not exist => no prior stats to load
    except IOError:
        pass

    # PROCESSING
    job_start = time()
    # process each job
    for job in config.jobs:
        print(f'Version {job.version}, {job.species}:')
        if job.run.skip:
            print(f'\tSKIPPED: {job.run.skip}\n')
            continue

        start = time()
        # build the database tables for the current job
        builder = DataBaseBuilder(job)                                          # type: DataBaseBuilder
        # store some information about the job if it successfully finished
        if not builder.cfg.run.skip:
            # add all stats provided by the current job in a way that it overrides similar information but not other
            # information that may be stored for this job (i.e. the original creation date)
            job_stats = builder.get_stats()
            for key, val in job_stats.items():
                stats[job.species][job.version][key] = val
            # delete old information fields no longer included in the stat output of the database to ensure consistent
            # formatting (except for the original creation date)
            for key in set(stats[job.species][job.version].keys()).difference(job_stats.keys()):
                if key != 'original date':
                    del stats[job.species][job.version][key]
            # if there is no information about the original creation date, assume that the current date is the first
            # time generating this particular job
            if 'original date' not in stats[job.species][job.version]:
                stats[job.species][job.version]['original date'] = stats[job.species][job.version]['update date']
        print(f'--> total: {time_stamp(start)}\n')

    # STAT OUTPUT
    with open(config.report_file, 'w') as file:
        # convert the defaultdict objects into normal dictionaries so that yaml can output them without error
        for key in set(stats.keys()):
            stats[key] = dict(stats[key])
        stats = dict(stats)
        # write the stat information
        yaml.safe_dump(stats, file)

    print(f'TOTAL: {time_stamp(job_start)}')


if __name__ == '__main__':
    generate()
