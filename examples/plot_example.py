from Fractional_Abundace.FractionalAbundance import FractionalAbundance
import matplotlib.pyplot as plt
from pathlib import Path

if __name__ == "__main__":
    FA = FractionalAbundance(element="Ar")
    FA.plot_FA_all()
    plt.savefig(Path(__file__).parents[1] / "images" / "FA_Ar_plot.png")
