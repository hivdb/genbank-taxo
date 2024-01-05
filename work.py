import requests
import re
from Bio import bgzf
import gzip
from pathlib import Path
import subprocess
import yaml
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


def check_organism(genbank_file, include_list=[], exclude_list=[]):

    # exclude has higher priority than include

    for i in exclude_list:
        if i.lower() in genbank_file.lower():
            return False

    if not include_list:
        return True

    for i in include_list:
        if i.lower() in genbank_file.lower():
            return True

    return False


def get_genbank_file_by_organism(gz_fd, include_list=[], exclude_list=[]):

    all_genbank_files = []
    one_genbank_file = []

    for line in gz_fd:
        if line.startswith('//'):
            one_genbank_file.append(line)
            all_genbank_files.append(''.join(one_genbank_file))
            one_genbank_file = []
        else:
            one_genbank_file.append(line)

    if one_genbank_file:
        print('show be no one genbank file')
        all_genbank_files.append(''.join(one_genbank_file))

    sel_genbank_files = [
        i
        for i in all_genbank_files
        if check_organism(i, include_list, exclude_list)
    ]

    return {
        'genbank_files': sel_genbank_files,
        'total': len(all_genbank_files),
    }


def download_genbank_files(ctx):
    gb_file_prefix = ctx['gb_file_prefix']
    release_path = ctx['release_path']

    file_name_list = get_file_name_list(gb_file_prefix)

    print('#Files', len(file_name_list))

    context = [
        (release_path, i)
        for i in file_name_list
    ]

    for i in concurrent_work(
            context, download_worker,
            multi_arg=True, progress=True):
        pass


def download_worker(release_path, gz_file_name):

    file_path = release_path / gz_file_name
    if file_path.exists():
        return True

    gz_file_url = GENBANK_FTP + gz_file_name
    download_gz(gz_file_url, release_path)

    return True


def select_genbank_files(ctx):
    release_path = ctx['release_path']

    include_list = ctx['include_list']
    exclude_list = ctx['exclude_list']

    dry_run = ctx['dry_run']

    save_path = ctx['save_path']
    save_path.mkdir(exist_ok=True, parents=True)

    total_seq = 0
    files = []

    for file_path in release_path.iterdir():
        if not file_path.name.endswith('.seq.gz'):
            continue

        files.append(file_path)

    files.sort(key=lambda x: int(re.search(r'\d+', x.name).group()))

    context = [
        (f, save_path, include_list, exclude_list, dry_run)
        for f in files
    ]

    for (num_seq, num_total) in concurrent_work(
            context, process_file, multi_arg=True, progress=True):
        # print('seq/total:', num_seq / num_total)
        total_seq += num_seq

    print(f'Selected #{total_seq} sequences')


def process_file(
        file_path, save_path, include_list, exclude_list, dry_run=True):
    bgzf_path = save_path / (file_path.stem + '.bgz')
    sel_path = save_path / (file_path.stem + '.sel.gz')

    with gzip.open(file_path, 'rt') as fd:
        seq_info = get_genbank_file_by_organism(
            fd, include_list, exclude_list)

        seq_list = seq_info['genbank_files']

    if not seq_list:
        return len(seq_list), seq_info['total']

    if dry_run:
        return len(seq_list), seq_info['total']

    with bgzf.BgzfWriter(bgzf_path, 'wb') as fd:
        for i in seq_list:
            fd.write(f'{i}\n')

    with bgzf.BgzfReader(bgzf_path, 'rb') as bgzfd:
        with gzip.open(sel_path, 'w') as gzfd:
            for i in bgzfd:
                gzfd.write(i)

    bgzf_path.unlink(missing_ok=True)

    return len(seq_list), seq_info['total']


def work():

    config = load_yaml('config.yml')
    config['db_path'] = Path(config['db_path']).expanduser().resolve()
    config['release_path'] = config['db_path'] / config['date']

    download_genbank_files(config)

    config['include_list'] = config.get('include_list', [])
    config['exclude_list'] = config.get('exclude_list', [])

    config['save_path'] = config['db_path'] / config['extract_folder']

    print('Include:', ', '.join(config['include_list']))
    print('Exclude:', ', '.join(config['exclude_list']))
    select_genbank_files(config)


if __name__ == '__main__':
    work()
