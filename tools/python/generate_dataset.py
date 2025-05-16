import random
from argparse import ArgumentParser
from math import ceil

def generate_target(alphabet, length):
    pattern = ''.join(random.choices(alphabet, k=length))

    return pattern

def generate_query(alphabet, pattern, error):
    nerrors = ceil(len(pattern) * error)

    query = pattern

    for i in range(nerrors):
        pos = random.randint(0, len(query) - 1)
        op = random.randint(0, 2)
        # Insertion
        if op == 0:
            query = query[:pos] + random.choice(alphabet) + query[pos:]
        # Deletion
        elif op == 1:
            query = query[:pos] + query[pos + 1:]
        # Substitution
        else:
            query = query[:pos] + \
                random.choice(alphabet.replace(query[pos], '')) + query[pos + 1:]

    return query


def generate_dataset(alphabet, npairs, lenght, error):
    for i in range(npairs):
        target = generate_target(alphabet, lenght)
        query = generate_query(alphabet, target, error)
        print('>', target, sep='')
        print('<', query, sep='')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--npairs', type=int, required=True,
                        help='Number of pairs to generate')
    parser.add_argument('-l', '--length', type=int, required=True,
                        help='Length of the target sequences')
    # parser.add_argument('-e', '--error', type=float, required=True,
    #                     help='Error rate')

    args = parser.parse_args()

    # if args.error < 0 or args.error > 1:
    #     raise ValueError('Error rate must be between 0 and 1')
    error = 0.2

    alphabet = 'ACGT'
    generate_dataset(alphabet, args.npairs, args.length, error)
