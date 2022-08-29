import requests
import re
from tqdm import tqdm
from Bio import SeqIO, bgzf
import gzip
from pathlib import Path
import subprocess
import yaml


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


def get_genbank_file_by_organism(gz_fd, organism_list=[]):

    genbank_file_list = []

    total = 0
    for rec in SeqIO.parse(gz_fd, 'genbank'):
        total += 1

        if rec.annotations['organism'].lower() in organism_list:
            genbank_file_list.append(rec)

    return {
        'genbank_files': genbank_file_list,
        'total': total,
    }


def select_genbank_files(ctx):
    gb_file_prefix = ctx['gb_file_prefix']
    base_path = ctx['base_path']
    organism_list = ctx['organism_list']

    file_name_list = get_file_name_list(gb_file_prefix)

    print('#Files', len(file_name_list))

    total_seq = 0

    for gz_file_name in tqdm(file_name_list):

        file_path = base_path / gz_file_name
        bgzf_path = base_path / (file_path.stem + '.bgz')
        sel_path = base_path / (file_path.stem + '.sel.gz')

        if sel_path.exists():
            with gzip.open(sel_path, 'rt') as fd:
                seq_list = list(SeqIO.parse(fd, 'genbank'))
                total_seq += len(seq_list)
                print(
                    'Get',
                    f"{len(seq_list)}",
                    f'total {total_seq}')
            continue

        file_path.unlink(missing_ok=True)

        gz_file_url = GENBANK_FTP + gz_file_name
        download_gz(gz_file_url, base_path)

        with gzip.open(file_path, 'rt') as fd:
            seq_info = get_genbank_file_by_organism(
                fd, organism_list)

            seq_list = seq_info['genbank_files']
            total_seq += len(seq_list)

            print(
                'Get',
                f"{len(seq_list)}/{seq_info['total']}",
                f'total {total_seq}')

        file_path.unlink(missing_ok=True)

        if not seq_list:
            continue

        with bgzf.BgzfWriter(bgzf_path, 'wb') as fd:
            SeqIO.write(seq_list, fd, 'genbank')

        with bgzf.BgzfReader(bgzf_path, 'rb') as bgzfd:
            with gzip.open(sel_path, 'w') as gzfd:
                for i in bgzfd:
                    gzfd.write(i)

        bgzf_path.unlink(missing_ok=True)


def work():

    config = load_yaml('config.yml')
    config['organism_list'] = [
        i.lower()
        for i in config['organism_list']
    ]
    config['base_path'] = Path(config['base_path']).expanduser().resolve()

    print('Select:', ', '.join(config['organism_list']))
    select_genbank_files(config)


if __name__ == '__main__':
    work()

# TODO check malform
