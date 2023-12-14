import logging
import os
import re

from src.definitions import DEFAULT_LOGGER_NAME, DEFAULT_LOGGER_LEVEL, DEFAULT_PTOOLS_CHAR_LIMIT
from src.lib.process import logging_helper


def read_delim_itr(fp, key_idx=0, val_indices=None, delim=None, skip=None):
    """Iterator to read delimited files
    Args:
        fp: Opened file object
        key_idx: The index for they key value
        val_indices: The indices for the values to retrieve
        delim: The list of substrings for delimiters
        skip: The list of substrings that will be skipped
    Raises:
    Yields: A key and it's values
     """
    if val_indices is None:
        val_indices = [1]
    if delim is None:
        delim = ["\t"]
    if skip is None:
        skip = ["!", "#"]
    delim_regex = '|'.join(map(re.escape, delim))
    for line in fp:
        if line.startswith(tuple(skip)):
            continue
        line = line.rstrip('\n')
        info = re.split(delim_regex, line)
        try:
            yield info[key_idx], [i for idx, i in enumerate(info) if idx in val_indices]
        except IndexError:
            continue


def read_e2p2_maps(ef_map_path, key_idx=0, val_idx=1, logging_level=DEFAULT_LOGGER_LEVEL,
                   logger_name=DEFAULT_LOGGER_NAME):
    """Read in mapping files of E2P2
    Args:
        ef_map_path: Path to an E2P2 mapping file
        key_idx: Index of the key
        val_idx: Index of the value
        logging_level: The logging level set for read map
        logger_name: The name of the logger for read map
    Raises:
    Returns:
        map_dict: Key to value of the E2P2 map
    """
    map_dict = {}
    logging_helper("Loading map: \"" + ef_map_path + "\"", logging_level=logging_level, logger_name=logger_name)
    with open(ef_map_path, 'r') as fp:
        for info in read_delim_itr(fp, key_idx=key_idx, val_indices=[val_idx]):
            if info:
                key = info[0]
                val = info[1]
                try:
                    map_dict[key] += val
                except KeyError:
                    map_dict.setdefault(key, val)
    return map_dict


def read_groups_by_start_itr(fp, start=None, skip=None):
    """Iterator to group file strings by their starting string
    Args:
        fp: Opened file object
        start: The list of substrings that indicates new groups
        skip: The list of substrings that will be skipped
    Raises:
    Yields: Group of lines based on 'start'
    """
    if start is None:
        start = ['>']
    if skip is None:
        skip = ["!", "#"]
    header, group = None, []
    for line in fp:
        line = line.rstrip('\n')
        if line.startswith(tuple(skip)):
            continue
        elif line.startswith(tuple(start)):
            if header:
                yield header, group
            header, group = line, []
        else:
            group.append(line)
    if header:
        yield header, group


def read_fasta(fp):
    """Iterator for reading fasta files. Source: Biopython
    Args:
        fp: file pointer to fasta file
    Raises:
    Yields:
        header: fasta header
        seq: fasta sequence
    """
    header, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith('>'):
            if header:
                yield header, '\n'.join(seq)
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield header, '\n'.join(seq)


def check_fasta_header(fasta_path, logger_name=DEFAULT_LOGGER_NAME):
    """Warn if fasta sequence ID length increases Pathway Tools current limit
    Args:
        fasta_path: The path to fasta input
        logger_name: The name of the logger for checking fasta header
    Raises: IndexError, KeyError
    Returns:
    """
    existing_headers = []
    with open(fasta_path, 'r') as fp:
        for header, seq in read_fasta(fp):
            try:
                header_info = re.split('[\s|]+', header)
                header_id = header_info[0].replace('>', '', 1)
                if len(header_id) > DEFAULT_PTOOLS_CHAR_LIMIT:
                    logging_helper("ID exceeds Pathway-Tools character limit: " + header_id,
                                   logging_level="WARNING", logger_name=logger_name)
                if header_id in existing_headers:
                    logging_helper("Duplicate IDs: " + header_id,
                                   logging_level="ERROR", logger_name=logger_name)
                    break
                else:
                    existing_headers.append(header_id)
            except (IndexError, KeyError):
                logging_helper("Cannot Parse Header: " + header,
                               logging_level="ERROR", logger_name=logger_name)
                continue


def remove_splice_variants_from_fasta(fasta_path, output_dir, prot_gene_map, logger_name=DEFAULT_LOGGER_NAME):
    """Remove splice variants from input fasta
    Args:
        fasta_path: Path to fasta input
        output_dir: Path to splice variants removed fasta output
        prot_gene_map: Path to Mapping file of protein IDs to gene IDs.
        logger_name: The name of the logger for remove splice variants
    Raises: IndexError, KeyError
    Returns:
    """
    logger = logging.getLogger(logger_name)
    fasta_dict = {}
    # The input's header should already be formatted
    file_name, file_extension = os.path.splitext(os.path.basename(fasta_path))
    prot_gene_map_dict = read_e2p2_maps(prot_gene_map, 0, 1)
    output_path = os.path.join(output_dir, file_name + '.rmspl' + file_extension)
    if os.path.isfile(output_path):
        logger.log(logging.WARNING, "Output path %s exists, will overwrite..." % output_path)
    with open(fasta_path, 'r') as fp:
        for header, seq in read_fasta(fp):
            try:
                header_info = re.split('[\s|]+', header)
                header_id = header_info[0].replace('>', '', 1)
                try:
                    locus = prot_gene_map_dict[header_id][0]
                except KeyError:
                    locus = header_id
                fasta_tuple = fasta_dict.setdefault(locus, (header, seq))
                if len(seq) > len(fasta_tuple[1]):
                    fasta_dict[locus] = (header, seq)
            except (IndexError, KeyError):
                logging_helper("Cannot Parse Header: " + header,
                               logging_level="WARNING", logger_name=logger_name)
                continue
    logging_helper("Removing splice variants from: \"" + fasta_path + "\"",
                   logging_level="INFO", logger_name=logger_name)
    with open(output_path, 'w') as op:
        for locus in sorted(fasta_dict.keys()):
            try:
                header = fasta_dict[locus][0]
                seq = fasta_dict[locus][1]
                op.write(header + '\n' + seq + '\n')
            except (IndexError, KeyError):
                continue
    return output_path


def get_all_seq_ids_from_fasta(fasta_path, logger_name=DEFAULT_LOGGER_NAME):
    """Get all sequence IDs from a fasta file
    Args:
        fasta_path: Path to fasta input
        logger_name: The name of the logger
    Raises: IndexError, KeyError
    Returns:
    """
    seq_ids = []
    with open(fasta_path, 'r') as fp:
        for header, seq in read_fasta(fp):
            try:
                header_info = re.split('[|\s]+', header)
                header_id = header_info[0].replace('>', '', 1)
                seq_ids.append(header_id)
            except (IndexError, KeyError):
                logging_helper("Cannot Parse Header: " + header,
                               logging_level="WARNING", logger_name=logger_name)
                continue
    return seq_ids
