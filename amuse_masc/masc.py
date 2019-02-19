from __future__ import print_function

try:
    from qtpy.uic import loadUi
    from qtpy.QtWidgets import QApplication, QWidget
    gui = "qt"
except ImportError:
    gui = "none"

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from amuse.lab import units, SeBa, write_set_to_file
from make_a_star_cluster import make_a_star_cluster


class MakeAClusterWidget(QWidget):

    def __init__(self):

        super(MakeAClusterWidget, self).__init__()

        self.length_unit = units.parsec

        # Load file made using Qt Designer
        self.widget = loadUi('masc.ui')
        plot_area = self.widget.plot_area

        self.fig = Figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.setup_axes()

        self.canvas = FigureCanvasQTAgg(self.fig)
        plot_area.addWidget(self.canvas)

    def setup_axes(self):
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
        self.ax.set_xlabel(self.length_unit)
        self.ax.set_ylabel(self.length_unit)
        self.ax.set_aspect(1)

    def create_cluster(self):
        number_of_stars = self.widget.number_of_stars_box.value()
        star_metallicity = self.widget.star_metallicity_box.value()
        mass_distribution = self.widget.mass_distribution_box.itemText(
                self.widget.mass_distribution_box.currentIndex(),
                )
        star_distribution = self.widget.star_distribution_box.itemText(
                self.widget.star_distribution_box.currentIndex(),
                )

        self.particles = make_a_star_cluster(
                initial_mass_function=mass_distribution,
                number_of_stars=number_of_stars,
                star_metallicity=star_metallicity,
                star_distribution=star_distribution,
                )

        self.evolve_cluster()
        self.widget.actionSave.setEnabled(True)

    def evolve_cluster(self):
        age = self.widget.age_box.value() | units.Myr
        SE = SeBa()
        SE.particles.add_particles(self.particles)
        SE.evolve_model(age)
        self.particles.radius = SE.particles.radius
        self.particles.luminosity = SE.particles.luminosity
        SE.stop()
        self.plot_cluster()

    def plot_cluster(self):
        x = self.particles.x.value_in(self.length_unit)
        y = self.particles.y.value_in(self.length_unit)

        s = 0.5 * self.particles.luminosity.value_in(units.LSun)**0.5

        c = "black"

        self.ax.cla()
        self.setup_axes()

        self.ax.scatter(x, y, s=s, facecolor=c, lw=0, alpha=0.5)
        self.canvas.draw()

    def save_file(self):
        clusterfile = "cluster.hdf5"
        write_set_to_file(self.particles, clusterfile, "amuse")
        print("Saved to %s" % clusterfile)


# Initialize application
app = QApplication([])


widget = MakeAClusterWidget()

# Connect events to buttons
widget.widget.createButton.clicked.connect(widget.create_cluster)
widget.widget.evolveButton.clicked.connect(widget.evolve_cluster)
widget.widget.actionSave.triggered.connect(widget.save_file)

widget.widget.show()

# Start 'event loop'
app.exec_()

widget = loadUi('masc.ui')
