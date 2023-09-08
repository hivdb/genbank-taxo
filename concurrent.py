from multiprocessing import Pool
from unittest import result  # noqa, TODO
from tqdm import tqdm
from more_itertools import ichunked


def concurrent_work(
        context_list, func,
        multi_arg=False, iter=False, POOL_SIZE=10,
        progress=False):

    worker = pool_map_work

    if iter:
        worker = iterative_work
    elif multi_arg:
        worker = pool_starmap_work

    if progress:
        return _with_bar(worker, context_list, func, multi_arg, POOL_SIZE)
    else:
        return _without_bar(worker, context_list, func, multi_arg, POOL_SIZE)


def _with_bar(worker, context_list, func, multi_arg, POOL_SIZE):

    with tqdm(total=len(context_list)) as pbar:
        for idx, ret in enumerate(worker(
                context_list, func, multi_arg=multi_arg, POOL_SIZE=POOL_SIZE)):
            idx = idx + 1
            if idx % POOL_SIZE == 0:
                pbar.update(POOL_SIZE)

            if idx == len(context_list):
                pbar.update(len(context_list) % POOL_SIZE)

            yield ret


def _without_bar(worker, context_list, func, multi_arg, POOL_SIZE):
    results = []
    for ret in worker(
            context_list, func, multi_arg=multi_arg, POOL_SIZE=POOL_SIZE):

        results.append(ret)

    return results


def iterative_work(context_list, func, **kwargs):

    multi_arg = kwargs['multi_arg']

    for c in context_list:
        if multi_arg:
            yield func(*c)
        else:
            yield func(c)


def pool_map_work(context_list, func, **kwargs):

    POOL_SIZE = kwargs['POOL_SIZE']

    with Pool(POOL_SIZE) as p:
        for ret in p.imap(func, context_list, POOL_SIZE):
            yield ret


def pool_starmap_work(context_list, func, **kwargs):

    POOL_SIZE = kwargs['POOL_SIZE']

    context_list = ichunked(context_list, POOL_SIZE)

    with Pool(POOL_SIZE) as p:
        for c in context_list:
            for ret in p.starmap(func, c):
                yield ret
