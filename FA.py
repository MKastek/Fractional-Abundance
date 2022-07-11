import argparse
from FractionalAbundance import FractionalAbundance


def init_parser():
    parser = argparse.ArgumentParser(description="Fractional Abundances")
    parser.add_argument("atom", help="name of the atom")
    parser.add_argument("SCD", help="name of the SCD coefficient file from Open ADAS", type=str)
    parser.add_argument("ACD", help="name of the ACD coefficient file from Open ADAS", type=str)
    return parser.parse_args()


if __name__ == '__main__':
    args = init_parser()
    FA = FractionalAbundance(atom=args.atom)
    FA.calculate()
    FA.plot_FA_all(i_Ne=40)