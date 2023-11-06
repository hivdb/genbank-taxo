import requests
import re
from tqdm import tqdm
from Bio import SeqIO, bgzf
import gzip
from pathlib import Path
import subprocess
import yaml
from datetime import datetime
from concurrent import concurrent_work


GENBANK_FTP = 'https://ftp.ncbi.nih.gov/genbank/'
RE_GB_GZ_FILE_NAME = r'([a-zA-Z]+)(\d+)\.seq\.gz'


def download_gz(url, folder):
    subprocess.run(["wget", url, '-P', folder, '-q'])


def load_yaml(file_path):
    return yaml.load(open(file_path), Loader=yaml.Loader)


def get_file_name_list(prefix):
    page = requests.get(GENBANK_FTP)
    file_names = list(set(re.findall(RE_GB_GZ_FILE_NAME, page.text)))

    file_names = [
        i
        for i in file_names
        if i[0] == prefix
    ]

    file_names.sort(key=lambda x: int(x[-1]))
    file_names = [
        f'{i[0]}{i[1]}.seq.gz'
        for i in file_names
    ]

    file_names = [
        i
        for i in file_names
        if i.startswith(prefix)
    ]

    return file_names


def get_genbank_file_by_organism(gz_fd, include_list=[], exclude_list=[]):

    genbank_file_list = []

    total = 0
    for rec in SeqIO.parse(gz_fd, 'genbank'):
        total += 1

        if include_list:
            if rec.annotations['organism'].lower() in include_list:
                genbank_file_list.append(rec)
            continue

        if exclude_list:
            if rec.annotations['organism'].lower() in exclude_list:
                continue

        genbank_file_list.append(rec)

    return {
        'genbank_files': genbank_file_list,
        'total': total,
    }


def download_genbank_files(ctx):
    gb_file_prefix = ctx['gb_file_prefix']
    base_path = ctx['base_path']

    file_name_list = get_file_name_list(gb_file_prefix)

    print('#Files', len(file_name_list))

    context = [
        (base_path, i)
        for i in file_name_list
    ]

    for i in concurrent_work(
            context, download_worker,
            multi_arg=True, progress=True):
        pass


def download_worker(base_path, gz_file_name):

    file_path = base_path / gz_file_name
    if file_path.exists():
        return

    gz_file_url = GENBANK_FTP + gz_file_name
    download_gz(gz_file_url, base_path)


def select_genbank_files(ctx):
    base_path = ctx['base_path']
    include_list = ctx['include_list']
    exclude_list = ctx['exclude_list']
    save_path = ctx['save_path']
    save_path.mkdir(exist_ok=True, parents=True)

    total_seq = 0
    files = []

    for file_path in base_path.iterdir():
        if not file_path.name.endswith('.seq.gz'):
            continue

        files.append(file_path)

    context = [
        (f, save_path, include_list, exclude_list)
        for f in files
    ]

    for i in concurrent_work(
            context, process_file, multi_arg=True, progress=True):
        print(i)


def process_file(file_path, save_path, include_list, exclude_list):
    bgzf_path = save_path / (file_path.stem + '.bgz')
    sel_path = save_path / (file_path.stem + '.sel.gz')

    with gzip.open(file_path, 'rt') as fd:
        seq_info = get_genbank_file_by_organism(
            fd, include_list, exclude_list)

        seq_list = seq_info['genbank_files']

    if not seq_list:
        return len(seq_list), seq_info['total']

    with bgzf.BgzfWriter(bgzf_path, 'wb') as fd:
        SeqIO.write(seq_list, fd, 'genbank')

    # TODO, why concurrent work not print progress bar?

    with bgzf.BgzfReader(bgzf_path, 'rb') as bgzfd:
        with gzip.open(sel_path, 'w') as gzfd:
            for i in bgzfd:
                gzfd.write(i)

    bgzf_path.unlink(missing_ok=True)

    return len(seq_list), seq_info['total']


def work():

    config = load_yaml('config.yml')
    config['include_list'] = [
        i.lower()
        for i in config['include_list']
    ]
    config['exclude_list'] = [
        i.lower()
        for i in config['exclude_list']
    ]
    config['db_path'] = Path(config['db_path']).expanduser().resolve()
    config['base_path'] = config['db_path'] / config['date']

    download_genbank_files(config)

    # config['base_path'] = (
    #     config['folder_path'] / (datetime.today().isoformat()[:10]))

    # config['save_path'] = config['db_path'] / 'hiv'

    # print('Select:', ', '.join(config['include_list']))
    # print('Exclude:', ', '.join(config['exclude_list']))
    # select_genbank_files(config)


if __name__ == '__main__':
    work()

# TODO check malform
