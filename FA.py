import argparse
from Fractional_Abundace.FractionalAbundance import FractionalAbundance


def init_parser():
    parser = argparse.ArgumentParser(description="Fractional Abundances")
    parser.add_argument("element", help="name of the atom")
    parser.add_argument(
        "SCD", help="name of the SCD coefficient file from Open ADAS", type=str
    )
    parser.add_argument(
        "ACD", help="name of the ACD coefficient file from Open ADAS", type=str
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = init_parser()
    FA = FractionalAbundance(element=args.element)
    FA.calculate()
    FA.plot_FA_all(index_Ne=40)
